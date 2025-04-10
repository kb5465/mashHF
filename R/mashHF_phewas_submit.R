#!/bin/bash

# Step 1: QT
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
  
# Step 1: BT
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
  
#######################
## SUBMIT STEP 2  QT ##
#######################

# Array of c(CHR, CHUNK_MAX) pairs
pairs=("1,2" "2,2" "3,2" "4,1" "5,2" "6,2" "7,2" "8,2" "9,2" "10,2" "11,2" "12,2" "13,2" "14,2" "15,2" "16,2" "17,2" "18,2" "19,2" "20,2" "21,2" "22,2")
pairs=("1,2")

# Loop through the pairs and execute the code for each pair
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


#######################
## SUBMIT STEP 2  BT ##
#######################

# Loop through the pairs and execute the code for each pair
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


#------------------------------#
#-------------TEST QT-------------#
#------------------------------#
CHR=1
CHUNK=1
CHUNK_MAX=197
dx run dxjupyterlab \
-iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_QTphewas_REGENIE_step2.sh" \
-icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
-ifeature="PYTHON_R" \
--name mashHF_QTphewas_step2_TEST \
--instance-type mem2_ssd1_v2_x4 \
--priority low \
--yes \
--destination exome-seq:/kbiddinger/projects/mashHF/results/step2/

#------------------------------#
#------------------------------#
  
  
#------------------------------#
#-------------TEST BT-------------#
#------------------------------#
CHR=1
CHUNK=1
CHUNK_MAX=2
dx run dxjupyterlab \
-iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_QTphewas_REGENIE_step2.sh" \
-icmd="bash *.sh ${CHR} ${CHUNK} ${CHUNK_MAX}" \
-ifeature="PYTHON_R" \
--name mashHF_BTphewas_step2_TEST \
--instance-type mem2_ssd1_v2_x4 \
--priority low \
--yes \
--destination exome-seq:/kbiddinger/projects/mashHF/results/step2/

#------------------------------#
#------------------------------#

  
  
dx run dxjupyterlab \
-iin="exome-seq:/kbiddinger/projects/mashHF/scripts/mashHF_extractBETASE_phewas.sh" \
-icmd="bash *.sh" \
-ifeature="PYTHON_R" \
--name mashHF_extractBETASE_phewas \
--instance-type mem2_ssd1_v2_x4 \
--priority low \
--yes \
--destination exome-seq:/kbiddinger/projects/mashHF/results/st_phewas

