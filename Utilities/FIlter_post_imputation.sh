#!/bin/bash
# This script processes imputed chunks (if you run chunking for imputation).
# This script will concatenate (merges 22 chromosome into a single vcf) and then filter the imputed variants based on the desired level of INFO score confidence.
# The end result is a sorted, filtered, and combined imputed VCF containing 22 chromosomes.
# this script requires BCFTOOLS to be installed.

mkdir Allchunk

for i in {1..22}; do 
bcftools concat chr$i/*.vcf.gz -Oz -o Allchunk/chr$i.allchunk.vcf.gz

# change the INFO threshold according to your needed level.
bcftools filter Allchunk/chr$i.allchunk.vcf.gz -i 'INFO>=0.8' -Oz -o Allchunk/chr$i.filtered.vcf.gz

bcftools sort Allchunk/chr$i.filtered.vcf.gz -Oz -o Allchunk/chr$i.filtered.sorted.vcf.gz

bcftools index Allchunk/chr$i.filtered.sorted.vcf.gz
done

bcftools concat Allchunk/*.filtered.vcf.gz -Oz -o Allchr.filtered.vcf.gz
bcftools index Allchr.filtered.vcf.gz

zcat Allchunk/*.allchunk.vcf.gz | grep -v "#" | wc -l
zcat Allchr.filtered.vcf.gz | grep -v "#" | wc -l