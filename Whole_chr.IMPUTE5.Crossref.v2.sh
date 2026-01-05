#!/bin/bash
#Modify the paths as needed
#usage: for i in {1..22}; do ./Whole_chr.IMPUTE5.Crossref.v2.sh $i; done

CHR=$1

mkdir -p IMPUTED/Input_Crossref/chr${CHR}/

./impute5_1.1.5_static \
    --h ~/REFPANELS/OA_CrossRef.chr${CHR}.imp5 \
    --m MAP/chr${CHR}.b37.gmap.gz \
    --g Phased/Input.phased.chr${CHR}.vcf.gz \
    --r ${CHR} --thread 35 \
    --o IMPUTED/Input_Crossref/chr${CHR}/chr${CHR}.imputed.vcf.gz \
    --l IMPUTED/Input_Crossref/chr${CHR}/chr${CHR}.imputed.log
