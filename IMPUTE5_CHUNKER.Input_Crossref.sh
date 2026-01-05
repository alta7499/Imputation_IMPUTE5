#!/bin/bash
#CHUKER
#Chunks chromosomal segments into segments for impute5 to run. final output allows impute5 to iterate through the arguments without further parsings. use the final output to the impute5 script.

CHR=$1

mkdir -p IMPUTED/Input_Crossref/chr${CHR}/

./imp5Chunker_1.1.5_static \
    --h ../REFPANELS/OA_Crossref/OA_GASP_SG10K.chr${CHR}.vcf.gz \
    --g Phased/Input.phased.chr${CHR}.vcf.gz \
    --r ${CHR} \
    --o IMPUTED/Input_Crossref/chr${CHR}.IMP5_coordinates.txt \
    --l IMPUTED/Input_Crossref/chr${CHR}.chunk.log

cat IMPUTED/Input_Crossref/chr${CHR}.IMP5_coordinates.txt | awk '{print $1, $2, $3, $4}' > IMPUTED/Input_CrossRef/TARGET_CHUNKS.chr${CHR}.txt

