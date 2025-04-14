#!/bin/bash

###############################################
# mashHF: Submit REGENIE Jobs via DNAnexus
# Author: Kiran Biddinger
# Description: Launches Step 1 and Step 2 REGENIE
# jobs for QT and BT traits, including tests and 
# post-processing extraction.
###############################################

#------------------------------------------------#
# STEP 1: RUN QT TRAITS
#------------------------------------------------#
dx run dxjupyterlab \
  -iin="exome-seq:kbiddinger/projects/mashHF/scripts/mashHF_QT_REGENIE_step1.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  -iduration=500 \
  --name mashHF_QT_step1 \
  --instance-type mem2_ssd1_v2_x96 \
  --priority normal \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/step1/

#------------------------------------------------#
# STEP 1: RUN BT TRAITS
#------------------------------------------------#
dx run dxjupyterlab \
  -iin="exome-seq:kbiddinger/projects/mashHF/scripts/mashHF_BT_REGENIE_step1.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  -iduration=500 \
  --name mashHF_BT_step1 \
  --instance-type mem2_ssd1_v2_x96 \
  --priority normal \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/step1/

#------------------------------------------------#
# STEP 2: SUBMIT QT TRAITS ACROSS CHROMOSOME CHUNKS
#------------------------------------------------#
pairs=("1,197" "2,118" "3,105" "4,72" "5,85" "6,101" "7,88" "8,65" "9,75" "10,71" \
       "11,127" "12,101" "13,31" "14,60" "15,57" "16,81" "17,113" "18,26" "19,140" \
       "20,54" "21,21" "22,43")

for pair in "${pairs[@]}"; do
  IFS="," read -r CHR CHUNK_MAX <<< "$pair"
  for (( CHUNK=1; CHUNK<=CHUNK_MAX; CHUNK++ )); do
    dx run dxjupyterlab \
      -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_QT_REGENIE_step2.sh" \
      -icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
      -ifeature="PYTHON_R" \
      --name mashHF_QT_step2_${CHR}_${CHUNK} \
      --instance-type mem2_ssd1_v2_x4 \
      --priority low \
      --yes \
      --destination exome-seq:/kbiddinger/projects/mashHF/results/step2/
  done
done

#------------------------------------------------#
# STEP 2: SUBMIT BT TRAITS ACROSS CHROMOSOME CHUNKS
#------------------------------------------------#
for pair in "${pairs[@]}"; do
  IFS="," read -r CHR CHUNK_MAX <<< "$pair"
  for (( CHUNK=1; CHUNK<=CHUNK_MAX; CHUNK++ )); do
    dx run dxjupyterlab \
      -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_BT_REGENIE_step2.sh" \
      -icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
      -ifeature="PYTHON_R" \
      --name mashHF_BT_step2_${CHR}_${CHUNK} \
      --instance-type mem2_ssd1_v2_x4 \
      --priority low \
      --yes \
      --destination exome-seq:/kbiddinger/projects/mashHF/results/step2/
  done
done

#------------------------------------------------#
# EXTRACT BETA AND SE AFTER STEP 2 COMPLETION
#------------------------------------------------#
dx run dxjupyterlab \
  -iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_extractBETASE.sh" \
  -icmd="bash *.sh" \
  -ifeature="PYTHON_R" \
  --name mashHF_extractBETASE \
  --instance-type mem2_ssd1_v2_x4 \
  --priority low \
  --yes \
  --destination exome-seq:/kbiddinger/projects/mashHF/results/
