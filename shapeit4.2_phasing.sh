#!/bin/bash
##This Script performes haplotype phasing for genotyping dataset. the desired input for this is a vcf, separated into individual Chromosomes. 
CHR=$1

##First, bcftools normalize and align alleles to a reference fasta. performed for each chr.
bcftools norm -c sw -f ~/Alvin/Tools/human_g1k_v37.fasta Input/Input.chr${CHR}.vcf.gz -Oz -o Input/Input.chr${CHR}.normed.vcf.gz --write-index

#Normalized dataset are then phased using shapeit, performed on each chr. Uses a Reference Panel phasing (recommmended for below 50 samples.
./shapeit4.2 --input Input/Input.chr${CHR}.normed.vcf.gz  --region ${CHR} --map ~/Alvin/IMPUTE5/MAP/chr${CHR}.b37.gmap.gz --output Phased/Input.phased.chr${CHR}.vcf.gz --thread 20 -H ~/dataset/OA_GASP_SG10K.chr${CHR}.vcf.gz ## -H flag is reference panel in vcf format. helps to phase accurately, especially in small number of samples

bcftools index ../Phased/Input.phased.chr${CHR}.vcf.gz
