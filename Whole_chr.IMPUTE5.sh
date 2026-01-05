#!/bin/bash
# Modify the paths as needed
# usage: for i in {1..22}; do ./Whole_chr.IMPUTE5.sh $i; done
# Requires: MAP files, Reference Panel (imp5 format), phased VCF
# This script is a continuation of the previous Phasing script.
# This is recommended for small files < 1M snps. for larger data, using a chunking protocol may be desireable to improve efficiency.
CHR=$1

mkdir -p IMPUTED/Input_refname/chr${CHR}/

./impute5_1.1.5_static \
    --h ./REFS/refname.chr${CHR}.imp5 \
    --m MAP/chr${CHR}.b37.gmap.gz \
    --g Phased/Input.phased.chr${CHR}.vcf.gz \
    --r ${CHR} --thread 35 \
    --o IMPUTED/Input_refname/chr${CHR}/chr${CHR}.imputed.vcf.gz \
    --l IMPUTED/Input_refname/chr${CHR}/chr${CHR}.imputed.log
