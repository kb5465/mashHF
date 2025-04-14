#!/bin/bash

###############################################
# mashHF: REGENIE Step 2 (QT Traits) - Chunked
# Author: Kiran Biddinger
# Description: Runs REGENIE step 2 rare variant
# association tests across chunks of gene sets
# with Cauchy p-value combination.
###############################################

# ----------------------------
# 1. Read Input Parameters
# ----------------------------
CHR=$1         # Chromosome number
CHUNK=$2       # Chunk index
CHUNK_TOTAL=$3 # Total number of chunks

# ----------------------------
# 2. Set Directories & Files
# ----------------------------
INPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/step1/"
OUTPUT_DIR_RAW="exome-seq:/kbiddinger/projects/mashHF/results/step2/rawassociation/"
OUTPUT_DIR_CAUCHY="exome-seq:/kbiddinger/projects/mashHF/results/step2/cauchycombination/"
PHENO_COVAR_FILE="/mnt/project/kbiddinger/projects/mashHF/data/pheno/mashHF_pheno.tsv"

# Create output directories on DNAnexus
dx mkdir -p ${OUTPUT_DIR_RAW}
dx mkdir -p ${OUTPUT_DIR_CAUCHY}

# ----------------------------
# 3. Download Tools and Scripts
# ----------------------------
# Download REGENIE
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.8.gz_x86_64_Linux_mkl.zip

# Download helper R packages and scripts
dx download exome-seq:/schoi/rpackages/rpackages4_1_3.tar.gz
tar -xvf rpackages4_1_3.tar.gz
git clone --branch patch-1 https://github.com/seanjosephjurgens/TOPMed_AFib_pipeline.git
git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git

# Download REGENIE step 1 LOCO files
dx download ${INPUT_DIR}mashHF_QT_stage1_out*.loco

# ----------------------------
# 4. Define Traits and Tissues
# ----------------------------
TRAIT_LIST=(NTproBNP Troponin LVMi LVESVi LVEDVi LVEF LVConc meanWT strain_rad strain_long strain_circ)
TRAITS=$(IFS=, ; echo "${TRAIT_LIST[*]}")

TISSUE_LIST=(all canonical Heart_Left_Ventricle_pext0.9 Muscle_Skeletal_pext0.9)
TISSUES=$(IFS=, ; echo "${TISSUE_LIST[*]}")

# ----------------------------
# 5. Filter Variants by Chunk and Tissue
# ----------------------------
chmod +x TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/filter_setinclusion_files.R
Rscript TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/filter_setinclusion_files.R \
  /mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.setinclusionfile.tsv \
  ${CHUNK} \
  ${CHUNK_TOTAL} \
  ${TISSUES} \
  UKBB_RVAT_groupingfile_chr${CHR}_chunk${CHUNK}.setinclusionfile.tsv

# ----------------------------
# 6. Run REGENIE Step 2 Association Testing
# ----------------------------
./regenie_v3.2.8.gz_x86_64_Linux_mkl \
  --step 2 \
  --pgen /mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c${CHR}_genotype_variant_sample_QCed \
  --covarFile ${PHENO_COVAR_FILE} \
  --covarColList age,age_squared,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,\
PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sex,batch \
  --catCovarList sex,batch \
  --phenoFile ${PHENO_COVAR_FILE} \
  --phenoColList ${TRAITS} \
  --qt \
  --pred /mnt/project/kbiddinger/projects/mashHF/results/step1/mashHF_QT_stage1_out_pred.list \
  --anno-file /mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.annotationfile.tsv \
  --set-list  /mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.setlistfile.tsv \
  --aaf-file  /mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.AAFfile.tsv \
  --extract-sets UKBB_RVAT_groupingfile_chr${CHR}_chunk${CHUNK}.setinclusionfile.tsv \
  --mask-def /mnt/project/kbiddinger/projects/LDL_RVAT/data/LDL_RVAT_groupingfile.maskdefinitionfile.tsv \
  --aaf-bins 0.01,0.001,0.00001 \
  --minMAC 5 \
  --vc-tests acato-full \
  --vc-maxAAF 0.01 \
  --vc-MACthr 10 \
  --joint sbat \
  --bsize 200 \
  --out mashHF_QT_stage2_out_rawassociationresults_chr${CHR}_${CHUNK}

# ----------------------------
# 7. Cauchy Combination for Each Trait
# ----------------------------
chmod +x TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/cauchy_construct_results.R

for TRAIT in "${TRAIT_LIST[@]}"; do
  regenie_outfile="mashHF_QT_stage2_out_rawassociationresults_chr${CHR}_${CHUNK}_${TRAIT}.regenie"
  cauchy_outfile="mashHF_QT_stage2_out_cauchycombinationresults_chr${CHR}_${CHUNK}_${TRAIT}.regenie"

  # Run Cauchy combination
  Rscript TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/cauchy_construct_results.R \
    ${regenie_outfile} \
    ${cauchy_outfile} \
    5 \
    0.01 \
    TRUE

  # Create directories for trait and chromosome
  dx mkdir -p ${OUTPUT_DIR_RAW}${TRAIT}/${CHR}/
  dx mkdir -p ${OUTPUT_DIR_CAUCHY}${TRAIT}/${CHR}/

  # Upload output files
  dx upload ${regenie_outfile} --destination ${OUTPUT_DIR_RAW}${TRAIT}/${CHR}/
  dx upload ${cauchy_outfile} --destination ${OUTPUT_DIR_CAUCHY}${TRAIT}/${CHR}/
done

# ----------------------------
# 8. Clean Up Working Directory
# ----------------------------
rm -rf /opt/notebooks/*
