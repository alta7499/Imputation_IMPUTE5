#!/bin/bash
# this script produces 2 collumn text file which contains SNP_ID and INFO score.
# the main purpose of this script is to generate an average of INFO score within an MAF bin. this tells us imputation quality across multiple MAF bins.
# RUN THIS SCRIPT INSIDE THE FOLDER YOU WANT TO RUN THE WHOLE PROCESS. NUMBER OF FILES WILL BE PRODUCED.

bfile=$1
echo "PlinkBinary: $bfile" > /dev/stderr


# make sure .impute2_info is within the folder.
for i in $(seq 1 22); do 
cat ../chr$i/*.impute2_info | awk '{ print "chr"'$i'":"$3"\t"$7 }' | grep -v "rs_id" > chr$i.infoscore
done

cat *.infoscore | grep -v "info" > ALLSNP_INFO_newID.txt

plink2 --bfile $bfile --set-all-var-ids chr@:# --make-bed --out $bfile.newID

plink2 --bfile $bfile.newID --max-maf 0.05 --make-bed --out 0_0.05.maf
cat 0_0.05.maf.bim | awk '{ print $2 }' > 0_0.05.maf.snplist.txt
range1=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0_0.05.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

plink2 --bfile $bfile.newID --maf 0.05 --max-maf 0.1 --make-bed --out 0.05_0.1.maf
cat 0.05_0.1.maf.bim |  awk '{ print $2 }' > 0.05_0.1.maf.snplist.txt
range2=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0.05_0.1.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

plink2 --bfile $bfile.newID --maf 0.1 --max-maf 0.2 --make-bed --out 0.1_0.2.maf
cat 0.1_0.2.maf.bim |  awk '{ print $2 }' > 0.1_0.2.maf.snplist.txt
range3=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0.1_0.2.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

plink2 --bfile $bfile.newID --maf 0.2 --max-maf 0.3 --make-bed --out 0.2_0.3.maf
cat 0.2_0.3.maf.bim |  awk '{ print $2 }' > 0.2_0.3.maf.snplist.txt
range4=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0.2_0.3.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

plink2 --bfile $bfile.newID --maf 0.3 --max-maf 0.4 --make-bed --out 0.3_0.4.maf
cat 0.3_0.4.maf.bim |  awk '{ print $2 }' > 0.3_0.4.maf.snplist.txt
range5=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0.3_0.4.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

plink2 --bfile $bfile.newID --maf 0.4 --max-maf 0.5 --make-bed --out 0.4_0.5.maf
cat 0.4_0.5.maf.bim |  awk '{ print $2 }' > 0.4_0.5.maf.snplist.txt
range6=$(cat ALLSNP_INFO_newID.txt | grep -Ff 0.4_0.5.maf.snplist.txt | awk '{ total += $2 } END { print total/NR }')

echo "Average bin score for maf 0-0.05: ${range1}"
echo "Average bin score for maf 0.05-0.1: ${range2}"
echo "Average bin score for maf 0.1-0.2: ${range3}"
echo "Average bin score for maf 0.2-0.3: ${range4}"
echo "Average bin score for maf 0.3-0.4: ${range5}"
echo "Average bin score for maf 0.4-0.5: ${range6}"
echo "${range1} ${range2} ${range3} ${range4} ${range5} ${range6}"