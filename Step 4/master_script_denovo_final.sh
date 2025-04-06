#!/bin/bash
# Script: master_script_denovo_final.sh
# Purpose: Coordinate annotation, filtering, DNM input generation, phasing prep, and phasing
# Usage: ./master_script_denovo_revised.sh <ACCOUNT> <FAMILY>

if [ $# -ne 2 ]; then
  echo "Error: Please provide ACCOUNT and FAMILY as arguments."
  echo "Usage: $0 <ACCOUNT> <FAMILY>"
  exit 1
fi

ACCOUNT=$1
FAMILY=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
SCRIPT_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/scripts"
PED_FILE="${BASE_DIR}/pedigree/${FAMILY}.ped"

[ -f "$PED_FILE" ] || { echo "Error: Pedigree file not found: $PED_FILE"; exit 1; }

SAMPLES=($(awk '{print $2}' "$PED_FILE"))
CHILDREN=($(awk '$3 != "0" && $4 != "0" {print $2}' "$PED_FILE"))
FATHER=$(awk '$3 == "0" && $5 == "1" {print $2}' "$PED_FILE" | head -n 1)
MOTHER=$(awk '$3 == "0" && $5 == "2" {print $2}' "$PED_FILE" | head -n 1)

[ -z "$FATHER" ] || [ -z "$MOTHER" ] && { echo "Error: Could not determine parents from $PED_FILE"; exit 1; }
[ ${#CHILDREN[@]} -eq 0 ] && { echo "Error: No children found in $PED_FILE"; exit 1; }

echo "Family: $FAMILY, Father: $FATHER, Mother: $MOTHER, Children: ${CHILDREN[*]}"
echo "Samples to process: ${SAMPLES[*]}"

# Step 1: Annotate trio with de novo mutations
echo "Running de novo annotation for $FAMILY..."
jid_anno=$(sbatch "${SCRIPT_DIR}/annotate_denovo_v2.sh" "$FAMILY" | awk '{print $4}')

# Step 2: Annotate known variants (trio)
jid_known=$(sbatch -A "$ACCOUNT" -J "${FAMILY}_annotate_known_v2" --dependency=afterok:$jid_anno "${SCRIPT_DIR}/annotate_known.sh" "$FAMILY" "trio" | awk '{print $4}')

# Step 3: Filter for de novo variants (per child)
declare -A filter_jobs
for CHILD in "${CHILDREN[@]}"; do
  jid=$(sbatch -A "$ACCOUNT" -J "${FAMILY}_filter_${CHILD}" --dependency=afterok:$jid_known "${SCRIPT_DIR}/filter_denovo_v3.sh" "$FAMILY" "$CHILD" | awk '{print $4}')
  filter_jobs[$CHILD]=$jid
done

# Step 4: Generate DNM input, and phase
for CHILD in "${CHILDREN[@]}"; do
  jid_dnm=$(sbatch -A "$ACCOUNT" -J "${FAMILY}_generate_dnm_${CHILD}" --dependency=afterok:${filter_jobs[$CHILD]} "${SCRIPT_DIR}/generate_dnm_input_v4.sh" "$FAMILY" "$CHILD" "$FATHER" "$MOTHER" | awk '{print $4}')
  jid_phase=$(sbatch -A "$ACCOUNT" -J "${FAMILY}_phase_denovo_${CHILD}" --dependency=afterok:$jid_dnm "${SCRIPT_DIR}/phase_denovo_v6.sh" "$FAMILY" "$CHILD" | awk '{print $4}')
done

echo "Monitor with: squeue -u rp793"
echo "Check logs in: ${BASE_DIR}/logs/"
