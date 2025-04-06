#!/bin/bash
# Script: phase_denovo_v6.sh
# Purpose: Run PhaseMyDeNovo to phase de novo variants (assumes phased CRAM)
# Usage: sbatch phase_denovo_v4.sh <FAMILY> <CHILD>

#SBATCH --ntasks=4
#SBATCH --mem=30000
#SBATCH --time=10:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake

FAMILY=$1
CHILD=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
#SBATCH -o ${BASE_DIR}/logs/phase_denovo_%j.out
#SBATCH -e ${BASE_DIR}/logs/phase_denovo_%j.err

if [ $# -ne 2 ]; then
  echo "Error: Please provide FAMILY and CHILD as arguments."
  echo "Usage: $0 <FAMILY> <CHILD>"
  exit 1
fi

PHASE_SCRIPT="/home/rp793/rds/hpc-work/tools/phase_my_denovos.py"
DNM_INPUT="${BASE_DIR}/vcf/annotated/${CHILD}_dnm_input.txt"
DNM_OUTPUT="${BASE_DIR}/vcf/annotated/${CHILD}_PhaseMyDeNovo_output.txt"

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load python/3.8

export PATH="$HOME/.local/bin:$PATH"  # Uses samtools 1.20 from manual install

[ -f "$PHASE_SCRIPT" ] || { echo "Error: PhaseMyDeNovo script not found: $PHASE_SCRIPT"; exit 1; }
[ -f "$DNM_INPUT" ] || { echo "Error: DNM input file not found: $DNM_INPUT"; exit 1; }

# Verify samtools is accessible
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found in PATH"; exit 1; }

echo "Phasing de novo variants for $CHILD..."
python3 "$PHASE_SCRIPT" -dnmfile "$DNM_INPUT" -outfile "$DNM_OUTPUT"

[ $? -eq 0 ] || { echo "Error: Phasing failed for $CHILD"; exit 1; }

echo "Phasing complete for $CHILD"
