#!/bin/bash

################################################################################
# Script: regenie_step1_phewas_QT.sh
# Author: Kiran Biddinger
# Project: mashHF
#
# Description:
#   REGENIE Step 1 for quantitative PHEWAS traits in UK Biobank data.
#
# Inputs:
#   - Genotype data: UKBB array data (bed/bim/fam)
#   - Phenotype + covariate file: mashHF_pheno_phewas.tsv
#
# Outputs:
#   - mashHF_QT_stage1_out* files (LOCO predictions, etc.)
#   - Uploaded to DNAnexus results/step1phewas directory
################################################################################

#------------------------------------------------------------------------------
# Create DNAnexus output directories
#------------------------------------------------------------------------------

dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/
dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/

#------------------------------------------------------------------------------
# Download and unzip REGENIE software
#------------------------------------------------------------------------------

wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.8.gz_x86_64_Linux_mkl.zip

#------------------------------------------------------------------------------
# Define traits and format for REGENIE
#------------------------------------------------------------------------------

TRAIT_LIST=(LDL_C BMI SBP DBP)
TRAITS=$(IFS=,; echo "${TRAIT_LIST[*]}")

#------------------------------------------------------------------------------
# Run REGENIE Step 1 (quantitative)
#------------------------------------------------------------------------------

./regenie_v3.2.8.gz_x86_64_Linux_mkl \
    --step 1 \
    --lowmem \
    --bed /mnt/project/genotyping_array/ukbb_genotyping_array_c1_22_merged \
    --extract /mnt/project/genotyping_array/array_snps_pass_strictMAF0.01.snplist \
    --covarFile /mnt/project/kbiddinger/projects/mashHF/data/pheno/phewas/mashHF_pheno_phewas.tsv \
    --covarColList age,age_squared,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sex,batch \
    --catCovarList sex,batch \
    --phenoFile /mnt/project/kbiddinger/projects/mashHF/data/pheno/phewas/mashHF_pheno_phewas.tsv \
    --phenoColList ${TRAITS} \
    --bsize 1000 \
    --qt \
    --out mashHF_QT_stage1_out

#------------------------------------------------------------------------------
# Upload results to DNAnexus and clean up
#------------------------------------------------------------------------------

dx upload mashHF_QT_stage1_out* --destination exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/
rm -rf *