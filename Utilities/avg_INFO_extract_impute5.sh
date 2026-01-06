#!/bin/bash
# this script produces 2 collumn text file which contains SNP_ID and INFO score.
# the main purpose of this script is to generate an average of INFO score within an MAF bin. this tells us imputation quality across multiple MAF bins.
# RUN THIS SCRIPT INSIDE THE FOLDER YOU WANT TO RUN THE WHOLE PROCESS. NUMBER OF FILES WILL BE PRODUCED.

file=$1
echo "VCF: $file" > /dev/stderr

# removes multiallelic from the imputed dataset
bcftools index $file.vcf.gz 
bcftools norm -m both $file.vcf.gz -Ou | bcftools view -M 2 -m 2 -Oz -o $file.no_MA.vcf.gz

# Removes monomorphic from the imputed dataset
bcftools view -i 'AF>0' $file.no_MA.vcf.gz -Oz -o $file.no_0_MAF.vcf.gz

INFO_0.05_below=$(zcat $file.no_0_MAF.vcf.gz | grep -v "#" | grep "IMP" | awk '{print $1, $2, $8}' | sed 's/=/\t/g' | sed 's/;/\t/g' | awk '$5<=0.05' | awk '{ total += $7 } END { print total/NR }')

INFO_0.05_above=$(zcat $file.no_0_MAF.vcf.gz | grep -v "#" | grep "IMP" | awk '{print $1, $2, $8}' | sed 's/=/\t/g' | sed 's/;/\t/g' | awk '$5>0.05 && $5<1' | awk '{ total += $7 } END { print total/NR }')

range1=$(bcftools view -i 'AF<=0.05' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range2=$(bcftools view -i 'AF>0.05 && AF <=0.1' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range3=$(bcftools view -i 'AF>0.1 && AF <=0.2' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range4=$(bcftools view -i 'AF>0.2 && AF <=0.3' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range5=$(bcftools view -i 'AF>0.3 && AF <=0.4' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range6=$(bcftools view -i 'AF>0.4 && AF <=0.5' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range7=$(bcftools view -i 'AF>0.5 && AF <1' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)
range8=$(bcftools view -i 'AF==1' $file.no_0_MAF.vcf.gz -Ou | grep -v "#" | wc -l)


NO_0=$(zcat $file.no_0_MAF.vcf.gz | grep -v "#" | wc -l)
noMA=$(zcat $file.no_MA.vcf.gz | grep -v "#" | wc -l)
TOTAL=$(zcat $file.vcf.gz | grep -v "#" | wc -l)
MA=$(($range3-$range2))
MAF_0=$(bcftools view -i 'AF==0' $file.vcf.gz -Ou | grep -v "#" | wc -l)
 
# print result and statistics
echo "AVG_INFO_MAF_0.05_below: $INFO_0.05_below"
echo "AVG_INFO_MAF_0.05_above: $INFO_0.05_above"
echo "MultiAllelic: $MA"
echo "NO_0: $NO_0"
echo "MAF_0: $MAF_0"

echo "Number of variants for maf 0-0.05: ${range1}"
echo "Number of variants for maf 0.05-0.1: ${range2}"
echo "Number of variants for maf 0.1-0.2: ${range3}"
echo "Number of variants for maf 0.2-0.3: ${range4}"
echo "Number of variants for maf 0.3-0.4: ${range5}"
echo "Number of variants for maf 0.4-0.5: ${range6}"
echo "Number of variants for maf above 0.5: ${range7}"
echo "Number of variants for homozygous Alternate: ${range8}"

echo "${range1} ${range2} ${range3} ${range4} ${range5} ${range6} ${range7} ${range8}"