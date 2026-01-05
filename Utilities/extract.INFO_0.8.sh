#!/bin/bash
#this script is specifically designed to list SNPS with certain INFO SCORE from IMPUTE2 results, which contain novel SNPs that requires to be renamed according to their chromosome number and position, instead of ".". 

for i in $(seq 1 22); do
cat ../../chr$i/*.impute2_info > ../../Allchunk/OA.chr$i.info
done

for i in $(seq 1 22); do 
cat ../../Allchunk/OA.chr$i.info | awk '$7 >= 0.8'| awk '{ print "chr"'$i'"\t"$2"\t"$3 }' | grep -v "rs_id" > ../../Allchunk/chr$i.pass.0.8.info
done
cat ../../Allchunk/*.pass.0.8.info > ../../Allchunk/pass.unprocessed.txt

grep "\." ../../Allchunk/pass.unprocessed.txt | awk '{ print $1"_"$3 }' > ../../Allchunk/novels.txt 
grep "rs" ../../Allchunk/pass.unprocessed.txt | awk '{ print $2 }' > ../../Allchunk/rsids.txt 
cat ../../Allchunk/novels.txt ../../Allchunk/rsids.txt > ../../Allchunk/ALL.0.8.snps.txt

../plink --bfile QCP.MAS.SSMPref --extract ../../Allchunk/ALL.0.8.snps.txt --make-bed --out filtered.QCP.MAS.SPMSref