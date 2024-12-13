#!/bin/bash

# Read input parameters
CHR=$1 
CHUNK=$2
CHUNK_TOTAL=$3

# Define input/output directories
INPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/step1/"
OUTPUT_DIR_RAW="exome-seq:/kbiddinger/projects/mashHF/results/step2/rawassociation/"
OUTPUT_DIR_CAUCHY="exome-seq:/kbiddinger/projects/mashHF/results/step2/cauchycombination/"

# Define phenotype file
PHENO_COVAR_FILE="/mnt/project/kbiddinger/projects/mashHF/data/pheno/mashHF_pheno.tsv"

# Make relevant directories
dx mkdir -p ${OUTPUT_DIR_RAW}
dx mkdir -p ${OUTPUT_DIR_CAUCHY}

# Set minimum allele count parameter
minMAC=5

# Download REGENIE software
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.8.gz_x86_64_Linux_mkl.zip

# Download helper scripts
dx download exome-seq:/schoi/rpackages/rpackages4_1_3.tar.gz
tar -xvf rpackages4_1_3.tar.gz
git clone --branch patch-1 https://github.com/seanjosephjurgens/TOPMed_AFib_pipeline.git
git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git

# Download results from step 1
dx download ${INPUT_DIR}mashHF_QT_stage1_out*.loco

# Define list of study traits
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

# Format list of traits
TRAITS=$(printf ",%s" "${TRAIT_LIST[@]}")
TRAITS=${TRAITS:1}

# Define tissues to keep (to create a tissue-specific subset of variants based on expression data)
TISSUE_LIST=(
    all
    canonical
    Heart_Left_Ventricle_pext0.9
    Muscle_Skeletal_pext0.9
)

# Format list of tissues
TISSUES=$(printf ",%s" "${TISSUE_LIST[@]}")
TISSUES=${TISSUES:1}

# Filter variants to current chromosome, chunk, and tissues of interest
chmod +x TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/filter_setinclusion_files.R
Rscript TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/filter_setinclusion_files.R  \
/mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.setinclusionfile.tsv  \
${CHUNK}  \
${CHUNK_TOTAL}  \
${TISSUES}  \
UKBB_RVAT_groupingfile_chr${CHR}_chunk${CHUNK}.setinclusionfile.tsv

# Run REGENIE step 2 of association testing
./regenie_v3.2.8.gz_x86_64_Linux_mkl \
    --step 2 \
    --pgen /mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c${CHR}_genotype_variant_sample_QCed \
    --covarFile ${PHENO_COVAR_FILE} \
    --covarColList age,age_squared,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sex,batch \
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
    --minMAC ${minMAC} \
    --vc-tests acato-full \
    --vc-maxAAF 0.01 \
    --vc-MACthr 10 \
    --joint sbat \
    --bsize 200 \
    --out mashHF_QT_stage2_out_rawassociationresults_chr${CHR}_${CHUNK}


# Make Cauchy combination helper script executable
chmod +x TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/cauchy_construct_results.R

# Iterate through traits
for TRAIT in ${TRAIT_LIST[@]}; do
    # Define output file names
    regenie_outfile=mashHF_QT_stage2_out_rawassociationresults_chr${CHR}_${CHUNK}_${TRAIT}.regenie
    cauchy_outfile=mashHF_QT_stage2_out_cauchycombinationresults_chr${CHR}_${CHUNK}_${TRAIT}.regenie

    # Call helper R script to perform Cauchy p-value combination
    Rscript TOPMed_AFib_pipeline/DNANexus/REGENIE_pipeline/v2/cauchy_construct_results.R  \
    ${regenie_outfile}  \
    ${cauchy_outfile}  \
    ${minMAC}  \
    0.01  \
    TRUE

    # Make output directories
    dx mkdir ${OUTPUT_DIR_RAW}${TRAIT}/
    dx mkdir ${OUTPUT_DIR_RAW}${TRAIT}/${CHR}/
    dx mkdir ${OUTPUT_DIR_CAUCHY}${TRAIT}/
    dx mkdir ${OUTPUT_DIR_CAUCHY}${TRAIT}/${CHR}/

    # Upload output
    dx upload ${regenie_outfile} --destination ${OUTPUT_DIR_RAW}${TRAIT}/${CHR}/
    dx upload ${cauchy_outfile} --destination ${OUTPUT_DIR_CAUCHY}${TRAIT}/${CHR}/
done

# Clean up workspace
rm -rf /opt/notebooks/*