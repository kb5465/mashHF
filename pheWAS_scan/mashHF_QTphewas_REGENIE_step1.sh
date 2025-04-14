#!/bin/bash

# Make relevant directories
dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/
dx mkdir -p exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/

# Download REGENIE software
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.8.gz_x86_64_Linux_mkl.zip

# Define list of study traits
TRAIT_LIST=(
    LDL_C
    BMI
    SBP
    DBP
)
TRAITS=$(printf ",%s" "${TRAIT_LIST[@]}")
TRAITS=${TRAITS:1}

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

dx upload mashHF_QT_stage1_out* --destination exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/

rm -rf *