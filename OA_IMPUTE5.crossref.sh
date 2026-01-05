#!/bin/bash
# This is a continuation of Chunker script. run IMPUTE5_CHUNKER.sh before running this.
# Requires: MAP files, Reference Panel (imp5 format), phased VCF, chunk files
# usage: for i in {1..22}; do cat TARGET_CHUNKS.chr$i.txt | while read Q; do ./Impute_Chunks_IMPUTE5.sh $Q done; done
# This is recommended for larger files > 1M snps. after running all the chunks to be imputed, make sure to check that there are no empty chunks. this may happen if the said chunk contains very less SNPs.
##In such case, you may want to increase chunk size, or altogether run a whole chromosome imputation. 

CHUNK=$1
CHR=$2
BUFFER=$3
TARGET=$4

./impute5_1.1.5_static \
    --h ./REFS/refname.chr${CHR}.imp5 \
    --m MAP/chr${CHR}.b37.gmap.gz \
    --g Phased/Input.phased.chr${CHR}.vcf.gz \
    --r ${TARGET} --buffer-region ${BUFFER} --thread 15 \
    --o IMPUTED/Input_refname/chr${CHR}/chr${CHR}.chunk${CHUNK}.imputed.vcf.gz \
    --l IMPUTED/Input_refname/chr${CHR}/chr${CHR}.chunk${CHUNK}.imputed.log
