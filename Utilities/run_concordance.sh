#!/usr/bin/env bash
# run_concordance.sh -- run `bcftools stats` on imputed vs truth VCF pairs and
# emit one .stats file per pair, ready for compile_concordance.py.
#
# By default both inputs are first reduced to biallelic SNPs
# (`bcftools view -m2 -M2 -v snps`); multiallelic sites are dropped, not split.
# See "Site preparation" in the usage text for how to change or skip that.
#
# Per-sample genotype concordance lives in the GCsS section, which bcftools only
# writes when (a) two files are given and (b) a samples option is supplied.
# so supplying sample is a MUST.
#
#   -s, --samples LIST       comma-separated list on the command line
#                            ("-" = all samples, "^A,B" = all except A and B)
#   -S, --samples-file FILE  the same list read from a file, one name per line
#
# Samples are matched between the two files BY NAME, so the IDs must be
# identical in the imputed and truth VCFs

set -euo pipefail
IMPUTED=""
TRUTH=""
MANIFEST=""
OUTDIR="bcftools_concordance"
LABEL="run"
SWAP=0
PREP=1
SITE_MODE="snps"
REF=""
CHECKREF="w"
RMDUP=0
PREPDIR=""
declare -a STATS_OPTS=()
declare -a REGION_OPTS=()
usage() {
    cat <<'EOF'
Usage:
  run_concordance.sh -i imputed.vcf.gz -t truth.vcf.gz [options]
  run_concordance.sh -m manifest.tsv [options]

Required (one of):
  -i FILE   imputed VCF/BCF (indexed)
  -t FILE   truth / "true-call" WGS VCF/BCF (indexed)
  -m FILE   manifest TSV: label <TAB> imputed <TAB> truth  (one pair per line,
            '#' comments allowed). Use this for per-chromosome or multi-run jobs.

Options:
  -o DIR    output directory                       [bcftools_concordance]
  -p STR    label for single-pair mode             [run]
  -r STR    regions (chr or chr:beg-end, comma-separated)  -> bcftools -r
  -R FILE   regions file (BED/tab)                          -> bcftools -R
  -T FILE   targets file, e.g. the masked sites only        -> bcftools -T
  -S FILE   samples file, one name per line (score a subset) -> bcftools -S
  -f EXPR   include expression, e.g. 'TYPE="snp"'           -> bcftools -i
  -e EXPR   exclude expression                              -> bcftools -e
  -k STR    --collapse mode (none|some|both|all|snps|indels) [none]
  -n INT    threads                                          [1]
  -X        swap file order (imputed first, truth second)
  -h        help

Site preparation (ON by default -- see notes):
  Default: both files are reduced to biallelic SNPs with
  `bcftools view -m2 -M2 -v snps`. Multiallelic records are DROPPED, not split.
  -P        skip preparation entirely, feed the raw files to stats
  -a        keep all variant types and multiallelics (no site filter)
  -M        split multiallelics with `norm -m -any` instead of dropping them
  -F FILE   reference FASTA (indexed); adds REF checking, and indel
            left-alignment when combined with -M
  -C MODE   --check-ref mode, needs -F: e|w|x|s                [w]
  -d        also drop duplicate records (--rm-dup exact, bcftools >= 1.12)
  -D DIR    where to write prepared files                      [OUTDIR/prep]

Notes:
  * Samples are matched by name, so the sample IDs must be identical in both
    VCFs. Only the intersection is scored; non-matching IDs drop out silently.
    Compare `bcftools query -l` on both files, and use `bcftools reheader -s`
    to rename if needed.
  * File order is truth-first by default. NRD and overall concordance are
    unaffected by the order, but the RR/RA/AA row labels are keyed off one file,
    so the per-genotype-class breakdown is not. Keep the order consistent.
  * Only sites present in BOTH files with matching alleles are compared. Sites
    the imputation dropped are counted as file-private, not as discordance.
  * Biallelic SNPs only, by default. Multiallelic sites are excluded from both
    files rather than split, because under -c none a record only matches one
    with an identical ALT string -- so a site the two callers represent
    differently would drop out of the comparison silently anyway. Excluding
    them up front makes that explicit and keeps the denominator honest.
    This also matches what GCsS scores: that section is SNPs only.
    Use -a to keep everything, -M to split instead of drop, -P to skip prep.
  * Prepared files are cached in OUTDIR/prep and reused across manifest rows.
EOF
}
require() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH" >&2; exit 1; }
}
check_pair() {
    local imputed="$1" truth="$2"
    for f in "$imputed" "$truth"; do
        [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
        if [[ ! -s "${f}.tbi" && ! -s "${f}.csi" ]]; then
            echo "WARN: no .tbi/.csi index for $f -- region options will fail" >&2
        fi
    done
    local n_shared n_truth n_imputed
    n_truth=$(bcftools query -l "$truth" | wc -l)
    n_imputed=$(bcftools query -l "$imputed" | wc -l)
    n_shared=$(comm -12 <(bcftools query -l "$truth" | sort) <(bcftools query -l "$imputed" | sort) | wc -l)
    if [[ "$n_shared" -eq 0 ]]; then
        echo "ERROR: no sample IDs shared between $truth and $imputed -- names must match exactly" >&2
        exit 1
    fi
    if [[ "$n_shared" -lt "$n_truth" || "$n_shared" -lt "$n_imputed" ]]; then
        echo "WARN: only ${n_shared} of ${n_truth} (truth) / ${n_imputed} (imputed) sample IDs match;" >&2
        echo "      the rest are excluded from concordance. Check for renamed IDs." >&2
    fi
    echo "  shared samples: ${n_shared}" >&2
}
prep_vcf() {
    # Prepares one input for comparison. Default (SITE_MODE=snps) keeps only
    # biallelic SNPs: `bcftools view -m2 -M2 -v snps`, which DROPS multiallelic
    # records outright rather than splitting them. Prints the output path on
    # stdout; progress goes to stderr. Cached by output name, so a truth VCF
    # shared by several manifest rows is prepared once.
    local src="$1" tag="$2"
    local base out
    base=$(basename "$src")
    base="${base%.gz}"
    base="${base%.vcf}"
    base="${base%.bcf}"
    base="${base}.${SITE_MODE}"
    if [[ "${#REGION_OPTS[@]}" -gt 0 ]]; then
        # different -r/-R subsets of the same source must not share a cache entry
        base="${base}.r$(printf '%s ' "${REGION_OPTS[@]}" | cksum | cut -d' ' -f1)"
    fi
    out="${PREPDIR}/${base}.bcf"
    if [[ -s "$out" && -s "${out}.csi" ]]; then
        echo "  ${tag}: reusing ${out}" >&2
        echo "$out"
        return
    fi
    echo "  ${tag}: preparing -> ${out}" >&2
    local -a view_opts=()
    local -a norm_opts=()
    if [[ "$SITE_MODE" == "snps" ]]; then view_opts+=(-m2 -M2 -v snps); fi
    if [[ "$SITE_MODE" == "split" ]]; then norm_opts+=(-m -any); fi
    if [[ -n "$REF" ]]; then norm_opts+=(-f "$REF" -c "$CHECKREF"); fi
    if [[ "$RMDUP" -eq 1 ]]; then norm_opts+=(--rm-dup exact); fi
    if [[ "${#norm_opts[@]}" -gt 0 ]]; then
        bcftools view -Ou ${REGION_OPTS[@]+"${REGION_OPTS[@]}"} ${view_opts[@]+"${view_opts[@]}"} "$src" \
            | bcftools norm "${norm_opts[@]}" --threads "$THREADS" -Ob -o "${out}.tmp"
    else
        bcftools view -Ob ${REGION_OPTS[@]+"${REGION_OPTS[@]}"} ${view_opts[@]+"${view_opts[@]}"} \
            --threads "$THREADS" -o "${out}.tmp" "$src"
    fi
    mv "${out}.tmp" "$out"
    bcftools index -f "$out"
    echo "$out"
}
run_one() {
    local label="$1" imputed="$2" truth="$3"
    local out="${OUTDIR}/${label}.stats"
    echo "[$(date +%T)] ${label}" >&2
    check_pair "$imputed" "$truth"
    if [[ "$PREP" -eq 1 ]]; then
        imputed=$(prep_vcf "$imputed" imputed)
        truth=$(prep_vcf "$truth" truth)
    fi
    local first="$truth" second="$imputed"
    if [[ "$SWAP" -eq 1 ]]; then first="$imputed"; second="$truth"; fi
    bcftools stats -s - "${STATS_OPTS[@]}" "$first" "$second" > "${out}.tmp"
    mv "${out}.tmp" "$out"
    if ! grep -q '^GCsS' "$out"; then
        echo "WARN: no GCsS records in ${out} -- no overlapping sites, or samples did not match" >&2
    fi
    echo "  wrote ${out}" >&2
}
COLLAPSE="none"
THREADS=1
while getopts "i:t:m:o:p:r:R:T:S:f:e:k:n:PaMNF:C:dD:Xh" opt; do
    case "$opt" in
        i) IMPUTED="$OPTARG" ;;
        t) TRUTH="$OPTARG" ;;
        m) MANIFEST="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        p) LABEL="$OPTARG" ;;
        r) REGION_OPTS+=(-r "$OPTARG") ;;
        R) REGION_OPTS+=(-R "$OPTARG") ;;
        T) STATS_OPTS+=(-T "$OPTARG") ;;
        S) STATS_OPTS+=(-S "$OPTARG") ;;
        f) STATS_OPTS+=(-i "$OPTARG") ;;
        e) STATS_OPTS+=(-e "$OPTARG") ;;
        k) COLLAPSE="$OPTARG" ;;
        n) THREADS="$OPTARG" ;;
        P) PREP=0 ;;
        a) SITE_MODE="all" ;;
        M) SITE_MODE="split" ;;
        N) : ;;
        F) REF="$OPTARG" ;;
        C) CHECKREF="$OPTARG" ;;
        d) RMDUP=1 ;;
        D) PREPDIR="$OPTARG" ;;
        X) SWAP=1 ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 1 ;;
    esac
done
STATS_OPTS+=(-c "$COLLAPSE" --threads "$THREADS")
STATS_OPTS+=(${REGION_OPTS[@]+"${REGION_OPTS[@]}"})
require bcftools
mkdir -p "$OUTDIR"
if [[ "$PREP" -eq 1 ]]; then
    PREPDIR="${PREPDIR:-${OUTDIR}/prep}"
    mkdir -p "$PREPDIR"
    echo "site mode: ${SITE_MODE}" >&2
    if [[ -n "$REF" ]]; then
        [[ -s "$REF" ]] || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
        [[ -s "${REF}.fai" ]] || echo "WARN: no ${REF}.fai -- run samtools faidx" >&2
    elif [[ "$SITE_MODE" == "split" ]]; then
        echo "WARN: -M without -F: multiallelics are split but indels are not" >&2
        echo "      left-aligned, so representation may still differ between files" >&2
    fi
fi
if [[ -n "$MANIFEST" ]]; then
    n=0
    while IFS=$'\t' read -r label imputed truth _rest; do
        [[ -z "${label:-}" || "${label:0:1}" == "#" ]] && continue
        run_one "$label" "$imputed" "$truth"
        n=$((n + 1))
    done < "$MANIFEST"
    echo "[$(date +%T)] done: ${n} pair(s)" >&2
elif [[ -n "$IMPUTED" && -n "$TRUTH" ]]; then
    run_one "$LABEL" "$IMPUTED" "$TRUTH"
else
    echo "ERROR: need -i and -t, or -m" >&2
    usage >&2
    exit 1
fi
echo "Next: compile_concordance.py ${OUTDIR}/*.stats -o ${OUTDIR}/concordance.tsv --aggregate" >&2
