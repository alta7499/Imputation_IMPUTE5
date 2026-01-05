#!/bin/bash
#RUN CHUNKER SCRIPT BEFORE THIS
#usage: for i in {1..22}; do cat TARGET_CHUNKS.chr$i.txt | while read Q; do ./OA_IMPUTE5.sh $Q done; done

CHUNK=$1
CHR=$2
BUFFER=$3
TARGET=$4

./impute5_1.1.5_static \
    --h ~/Alvin/REFPANELS/OA_Crossref.chr${CHR}.imp5 \
    --m MAP/chr${CHR}.b37.gmap.gz \
    --g Phased/Input.phased.chr${CHR}.vcf.gz \
    --r ${TARGET} --buffer-region ${BUFFER} --thread 15 \
    --o IMPUTED/Input_Crossref/chr${CHR}/chr${CHR}.chunk${CHUNK}.imputed.vcf.gz \
    --l IMPUTED/Input_Crossref/chr${CHR}/chr${CHR}.chunk${CHUNK}.imputed.log
