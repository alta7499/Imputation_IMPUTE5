#!/usr/bin/env python3
"""Compile per-sample genotype concordance from `bcftools stats` output.

Parses the GCsS ("genotype concordance by sample") section of one or more .stats
files and writes a tidy TSV: one row per (label, section, sample) with raw
match/mismatch counts plus derived rates.

Column positions are read from the '# GCsS [2]id [3]sample ...' header line in
each file rather than hard-coded, so this survives bcftools version differences.

Usage:
  compile_concordance.py out/*.stats -o concordance.tsv --aggregate
"""
import argparse
import os
import sys
# (output name, tokens that must appear, tokens that must NOT appear)
COUNT_SPECS = [
    ("rr_match", ("rr", "match"), ("mismatch",)),
    ("ra_match", ("ra", "match"), ("mismatch",)),
    ("aa_match", ("aa", "match"), ("mismatch",)),
    ("rr_mismatch", ("rr", "mismatch"), ()),
    ("ra_mismatch", ("ra", "mismatch"), ()),
    ("aa_mismatch", ("aa", "mismatch"), ()),
]
COUNT_NAMES = [spec[0] for spec in COUNT_SPECS]
OUT_COLS = ["label", "section", "sample", "n_compared"] + COUNT_NAMES + [
    "concordance", "nrd", "nrc", "homref_concordance", "het_concordance",
    "homalt_concordance", "reported_nrd", "dosage_r2",
]

def parse_header(line):
    """'# GCsS\t[2]id\t[3]sample\t...' -> (section, {normalised name: field index})."""
    fields = line.rstrip("\n").split("\t")
    section = fields[0].lstrip("#").strip()
    cols = {}
    for field in fields[1:]:
        field = field.strip()
        if not field.startswith("["):
            continue
        idx_str, _, name = field[1:].partition("]")
        try:
            cols[name.strip().lower()] = int(idx_str) - 1
        except ValueError:
            continue
    return section, cols

def find_col(cols, need, forbid):
    for name, idx in cols.items():
        if all(t in name for t in need) and not any(t in name for t in forbid):
            return idx
    return None

def build_layout(cols):
    """Map the GCsS header to field indices we care about; None if unusable."""
    layout = {"sample": cols.get("sample")}
    if layout["sample"] is None:
        return None
    for out_name, need, forbid in COUNT_SPECS:
        idx = find_col(cols, need, forbid)
        if idx is None:
            return None
        layout[out_name] = idx
    layout["reported_nrd"] = find_col(cols, ("discordance",), ())
    layout["dosage_r2"] = find_col(cols, ("dosage",), ())
    return layout

def parse_stats_file(path):
    """Yield dicts of raw counts for every GC*S-style per-sample section."""
    layouts = {}
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                if "\t[" not in line:
                    continue
                section, cols = parse_header(line)
                if not section.startswith("GC") or "sample" not in cols:
                    continue
                layout = build_layout(cols)
                if layout is not None:
                    layouts[section] = layout
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            section = fields[0]
            layout = layouts.get(section)
            if layout is None:
                continue
            row = {"section": section, "sample": fields[layout["sample"]]}
            try:
                for name in COUNT_NAMES:
                    row[name] = int(float(fields[layout[name]]))
            except (IndexError, ValueError):
                sys.stderr.write("WARN: skipping malformed %s line in %s\n" % (section, path))
                continue
            for name in ("reported_nrd", "dosage_r2"):
                idx = layout[name]
                row[name] = fields[idx] if idx is not None and idx < len(fields) else ""
            rows.append(row)
    return rows

def ratio(num, den):
    return num / den if den > 0 else float("nan")

def derive(row):
    """Add rates to a row of raw counts. Definitions:
      concordance = all matches / all compared genotypes
      NRD         = all mismatches / (compared genotypes - RR matches)
      NRC         = 1 - NRD
    NRD is symmetric in file order; the per-class rates below are not -- they are
    conditional on whichever file bcftools used for the row labels.
    """
    counts = {name: row[name] for name in COUNT_NAMES}
    matches = counts["rr_match"] + counts["ra_match"] + counts["aa_match"]
    mismatches = counts["rr_mismatch"] + counts["ra_mismatch"] + counts["aa_mismatch"]
    total = matches + mismatches
    row["n_compared"] = total
    row["concordance"] = ratio(matches, total)
    row["nrd"] = ratio(mismatches, total - counts["rr_match"])
    row["nrc"] = 1.0 - row["nrd"] if row["nrd"] == row["nrd"] else float("nan")
    row["homref_concordance"] = ratio(counts["rr_match"], counts["rr_match"] + counts["rr_mismatch"])
    row["het_concordance"] = ratio(counts["ra_match"], counts["ra_match"] + counts["ra_mismatch"])
    row["homalt_concordance"] = ratio(counts["aa_match"], counts["aa_match"] + counts["aa_mismatch"])
    return row

def aggregate(rows, label="ALL"):
    """Sum counts across labels, per (section, sample). Rates are recomputed from
    the summed counts -- never averaged. Reported NRD and dosage r2 are dropped
    because they cannot be pooled."""
    merged = {}
    for row in rows:
        key = (row["section"], row["sample"])
        target = merged.setdefault(key, {
            "label": label, "section": row["section"], "sample": row["sample"],
            "reported_nrd": "", "dosage_r2": "",
            **{name: 0 for name in COUNT_NAMES},
        })
        for name in COUNT_NAMES:
            target[name] += row[name]
    return [derive(row) for row in merged.values()]

def fmt(value):
    if isinstance(value, float):
        return "NA" if value != value else "%.6f" % value
    return str(value)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stats", nargs="+", help="one or more bcftools .stats files")
    parser.add_argument("-o", "--out", default="-", help="output TSV ('-' for stdout)")
    parser.add_argument("--aggregate", action="store_true",
                        help="also emit an 'ALL' row per sample, summing counts across inputs "
                             "(use when inputs are chromosomes/chunks of one run)")
    parser.add_argument("--aggregate-only", action="store_true", help="emit only the summed rows")
    parser.add_argument("--quiet", action="store_true", help="suppress the stderr summary")
    args = parser.parse_args()
    rows = []
    for path in args.stats:
        label = os.path.basename(path)
        for suffix in (".stats", ".txt", ".tsv"):
            if label.endswith(suffix):
                label = label[: -len(suffix)]
        parsed = parse_stats_file(path)
        if not parsed:
            sys.stderr.write("WARN: no per-sample concordance section in %s "
                             "(need two input VCFs and -s -)\n" % path)
        for row in parsed:
            row["label"] = label
            rows.append(derive(row))
    if not rows:
        sys.stderr.write("ERROR: nothing to compile\n")
        return 1
    out_rows = [] if args.aggregate_only else list(rows)
    if args.aggregate or args.aggregate_only:
        out_rows.extend(aggregate(rows))
    handle = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        handle.write("\t".join(OUT_COLS) + "\n")
        for row in out_rows:
            handle.write("\t".join(fmt(row.get(col, "")) for col in OUT_COLS) + "\n")
    finally:
        if handle is not sys.stdout:
            handle.close()
    if not args.quiet:
        finite = [r["nrd"] for r in out_rows if r["nrd"] == r["nrd"]]
        conc = [r["concordance"] for r in out_rows if r["concordance"] == r["concordance"]]
        sys.stderr.write("samples x runs: %d\n" % len(out_rows))
        if finite:
            finite.sort()
            sys.stderr.write("NRD          median %.5f  min %.5f  max %.5f\n"
                             % (finite[len(finite) // 2], finite[0], finite[-1]))
        if conc:
            conc.sort()
            sys.stderr.write("concordance  median %.5f  min %.5f  max %.5f\n"
                             % (conc[len(conc) // 2], conc[0], conc[-1]))
    return 0

if __name__ == "__main__":
    sys.exit(main())
