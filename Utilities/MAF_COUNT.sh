#!/bin/bash
# this script counts variants within each maf bins
# usage: ./MAF_COUNT.sh [plink_bin]

bfile=$1
echo "PlinkBinary: $bfile" > /dev/stderr

./plink --bfile $bfile --freq --out $bfile.MAF

range1=$(cat $bfile.MAF.frq | awk '$5 >= 0' | awk '$5 <= 0.05' | wc -l)
range2=$(cat $bfile.MAF.frq | awk '$5 > 0.05' | awk '$5 <= 0.1' | wc -l)
range3=$(cat $bfile.MAF.frq | awk '$5 > 0.1' | awk '$5 <= 0.2' | wc -l)
range4=$(cat $bfile.MAF.frq | awk '$5 > 0.2' | awk '$5 <= 0.3' | wc -l)
range5=$(cat $bfile.MAF.frq | awk '$5 > 0.3' | awk '$5 <= 0.4' | wc -l)
range6=$(cat $bfile.MAF.frq | awk '$5 > 0.4' | awk '$5 <= 0.5' | wc -l)

echo "Number of variants for maf 0-0.05: ${range1}"
echo "Number of variants for maf 0.05-0.1: ${range2}"
echo "Number of variants for maf 0.1-0.2: ${range3}"
echo "Number of variants for maf 0.2-0.3: ${range4}"
echo "Number of variants for maf 0.3-0.4: ${range5}"
echo "Number of variants for maf 0.4-0.5: ${range6}"
echo "${range1} ${range2} ${range3} ${range4} ${range5} ${range6}"
