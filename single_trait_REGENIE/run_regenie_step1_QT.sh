#!/bin/bash

###############################################
# mashHF: REGENIE Step 1 - QT Traits
# Author: Kiran Biddinger
# Description: Run REGENIE step 1 on quantitative
# traits using UK Biobank array data.
###############################################

# ----------------------------------
# 1. Create Output Directories on DNAnexus
# ----------------------------------
dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step1/
dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step2/rawassociation/

# ----------------------------------
# 2. Download and Unzip REGENIE Binary
# ----------------------------------
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.8.gz_x86_64_Linux_mkl.zip

# ----------------------------------
# 3. Define Traits for Analysis
# ----------------------------------
TRAIT_LIST=(
    NTproBNP
    Troponin
    LVMi
    LVESVi
    LVEDVi
    LVEF
    LVConc
    meanWT
    strain_rad
    strain_long
    strain_circ
)

# Convert array to comma-separated string
TRAITS=$(printf ",%s" "${TRAIT_LIST[@]}")
TRAITS=${TRAITS:1}  # Remove leading comma

# ----------------------------------
# 4. Run REGENIE Step 1
# ----------------------------------
./regenie_v3.2.8.gz_x86_64_Linux_mkl \
  --step 1 \
  --lowmem \
  --bed /mnt/project/genotyping_array/ukbb_genotyping_array_c1_22_merged \
  --extract /mnt/project/genotyping_array/array_snps_pass_strictMAF0.01.snplist \
  --covarFile /mnt/project/kbiddinger/projects/mashHF/data/pheno/mashHF_pheno.tsv \
  --covarColList age,age_squared,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,\
PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sex,batch \
  --catCovarList sex,batch \
  --phenoFile /mnt/project/kbiddinger/projects/mashHF/data/pheno/mashHF_pheno.tsv \
  --phenoColList ${TRAITS} \
  --bsize 1000 \
  --qt \
  --out mashHF_QT_stage1_out

# ----------------------------------
# 5. Upload Output Files to DNAnexus
# ----------------------------------
dx upload mashHF_QT_stage1_out* \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/step1/

# ----------------------------------
# 6. Clean Working Directory
# ----------------------------------
rm -rf *