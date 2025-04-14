#!/bin/bash

# Read input parameters
CHR=$1 
CHUNK=$2
CHUNK_TOTAL=$3

# Define input/output directories
INPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/"
OUTPUT_DIR_RAW="exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/"
OUTPUT_DIR_CAUCHY="exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/cauchycombination/"

# Define phenotype file
PHENO_COVAR_FILE="/mnt/project/kbiddinger/projects/mashHF/data/pheno/phewas/mashHF_pheno_phewas.tsv"

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
    LDL_C
    BMI
    SBP
    DBP
)

# Format list of traits
TRAITS=$(printf ",%s" "${TRAIT_LIST[@]}")
TRAITS=${TRAITS:1}

# Define tissues to keep (to create a tissue-specific subset of variants based on expression data)
TISSUE_LIST=(
    canonical
)

# Format list of tissues
TISSUES=$(printf ",%s" "${TISSUE_LIST[@]}")
TISSUES=${TISSUES:1}

# Print list of relevant genes
GENES_LIST="ENSG00000155657,ENSG00000134571,ENSG00000006747,ENSG00000007306,ENSG00000041515,ENSG00000044459,ENSG00000046889,ENSG00000048052,ENSG00000068781,ENSG00000069248,ENSG00000072210,ENSG00000080007,ENSG00000083635,ENSG00000086506,ENSG00000088827,ENSG00000099899,ENSG00000100344,ENSG00000100360,ENSG00000100441,ENSG00000101331,ENSG00000103150,ENSG00000105122,ENSG00000105136,ENSG00000105671,ENSG00000106328,ENSG00000108515,ENSG00000109736,ENSG00000109944,ENSG00000112280,ENSG00000115239,ENSG00000115525,ENSG00000115705,ENSG00000120279,ENSG00000123064,ENSG00000129696,ENSG00000132463,ENSG00000132635,ENSG00000132781,ENSG00000133048,ENSG00000133103,ENSG00000133454,ENSG00000133962,ENSG00000134014,ENSG00000134070,ENSG00000136273,ENSG00000136783,ENSG00000137198,ENSG00000138346,ENSG00000139988,ENSG00000140400,ENSG00000142082,ENSG00000142192,ENSG00000145949,ENSG00000146282,ENSG00000146350,ENSG00000148399,ENSG00000149054,ENSG00000151779,ENSG00000154102,ENSG00000157870,ENSG00000160191,ENSG00000161958,ENSG00000164404,ENSG00000164483,ENSG00000165071,ENSG00000165115,ENSG00000166529,ENSG00000167515,ENSG00000167733,ENSG00000169903,ENSG00000170379,ENSG00000170445,ENSG00000171564,ENSG00000176485,ENSG00000178665,ENSG00000179044,ENSG00000179142,ENSG00000179270,ENSG00000179709,ENSG00000181031,ENSG00000183844,ENSG00000188959,ENSG00000196188,ENSG00000198483,ENSG00000203780,ENSG00000205129,ENSG00000213563,ENSG00000214814,ENSG00000235718,ENSG00000242441,ENSG00000248385,ENSG00000248919,ENSG00000249222,ENSG00000253910,ENSG00000274391,ENSG00000278558,ENSG00000129151,ENSG00000215568,ENSG00000099957"

# Filter variants to current chromosome, chunk, and tissues of interest
chmod +x /mnt/project/kbiddinger/projects/mashHF/scripts/filter_setinclusion_files.R
Rscript /mnt/project/kbiddinger/projects/mashHF/scripts/filter_setinclusion_files.R  \
/mnt/project/exome_450k_plink/merged/annotation/group/REGENIE/tx_annotate/hclof_missense/all__canonical__tx_annotate_complete/UKBB_RVAT_groupingfile_chr${CHR}.setinclusionfile.tsv  \
${CHUNK}  \
${CHUNK_TOTAL}  \
${GENES_LIST} \
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
    --pred /mnt/project/kbiddinger/projects/mashHF/results/step1phewas/mashHF_QT_stage1_out_pred.list \
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