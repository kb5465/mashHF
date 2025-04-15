#!/bin/bash

################################################################################
# Script: submit_phewas_jobs.sh
# Author: Kiran Biddinger
# Project: mashHF
#
# Description:
#   Launches REGENIE step 1 and step 2 jobs for quantitative and binary traits
#   in the phewas extension of the mashHF pipeline. Also includes test jobs and
#   final extraction of effect estimates.
#
# Usage:
#   bash submit_phewas_jobs.sh
#
# Output:
#   Results will be written to:
#     exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/
#     exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/
#     exome-seq:/kbiddinger/projects/mashHF/results/st_phewas/
################################################################################

##################################
# Step 1: Quantitative Traits (QT)
##################################

dx run dxjupyterlab \
  -iin="exome-seq:kbiddinger/projects/mashHF/scripts/mashHF_QTphewas_REGENIE_step1.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  -iduration=500 \
  --name mashHF_QTphewas_step1 \
  --instance-type mem2_ssd1_v2_x96 \
  --priority normal \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/

#################################
# Step 1: Binary Traits (BT)
#################################

dx run dxjupyterlab \
  -iin="exome-seq:kbiddinger/projects/mashHF/scripts/mashHF_BTphewas_REGENIE_step1.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  -iduration=500 \
  --name mashHF_BTphewas_step1 \
  --instance-type mem2_ssd1_v2_x96 \
  --priority normal \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/step1phewas/

#########################
# Step 2: QT - Full Loop
#########################

# Array of chromosome/chunk pairs
pairs=("1,2" "2,2" "3,2" "4,1" "5,2" "6,2" "7,2" "8,2" "9,2" "10,2" "11,2" "12,2" \
       "13,2" "14,2" "15,2" "16,2" "17,2" "18,2" "19,2" "20,2" "21,2" "22,2")

for pair in "${pairs[@]}"; do
  IFS="," read -r CHR CHUNK_MAX <<< "$pair"
  for (( CHUNK=1; CHUNK<=CHUNK_MAX; CHUNK++ )); do
    dx run dxjupyterlab \
      -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_QTphewas_REGENIE_step2.sh" \
      -icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
      -ifeature="PYTHON_R" \
      --name mashHF_QTphewas_step2_${CHR}_${CHUNK} \
      --instance-type mem2_ssd1_v2_x4 \
      --priority low \
      --yes \
      --destination exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/
  done
done

#########################
# Step 2: BT - Full Loop
#########################

for pair in "${pairs[@]}"; do
  IFS="," read -r CHR CHUNK_MAX <<< "$pair"
  for (( CHUNK=1; CHUNK<=CHUNK_MAX; CHUNK++ )); do
    dx run dxjupyterlab \
      -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_BTphewas_REGENIE_step2.sh" \
      -icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
      -ifeature="PYTHON_R" \
      --name mashHF_BTphewas_step2_${CHR}_${CHUNK} \
      --instance-type mem2_ssd1_v2_x4 \
      --priority low \
      --yes \
      --destination exome-seq:/kbiddinger/projects/mashHF/results/step2phewas/
  done
done

##########################
# Final: Extract BETAs + SEs
##########################

dx run dxjupyterlab \
  -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_extractBETASE_phewas.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  --name mashHF_extractBETASE_phewas \
  --instance-type mem2_ssd1_v2_x4 \
  --priority low \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/st_phewas/