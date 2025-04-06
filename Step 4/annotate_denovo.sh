#!/bin/bash
# Script: annotate_denovo_v2.sh
# Purpose: Annotate a trio VCF with PossibleDeNovo using GATK VariantAnnotator
# Usage: sbatch annotate_denovo.sh <FAMILY> (e.g., sbatch annotate_denovo.sh GS-1)
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/%x/logs/%xdenovo_annotate%j.out
#SBATCH -e /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/%x/logs/%xdenovo_annotate%j.err

if [ $# -ne 1 ]; then
  echo "Error: Please provide FAMILY as argument."
  echo "Usage: $0 <FAMILY>"
  exit 1
fi

# Load required modules
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0
module load bcftools/1.19/gcc/3thngzmz

FAMILY=$1
ACCOUNT=${SLURM_JOB_ACCOUNT:-"CAGAN-SL2-CPU"}
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"
PED_FILE="${BASE_DIR}/pedigree/${FAMILY}.ped"

# Ensure output directories exist
mkdir -p "${BASE_DIR}/vcf/annotated" "${BASE_DIR}/logs"

# Define VCF and log paths (trio-based)
INPUT_VCF="${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz"
OUTPUT_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.vcf.gz"

# Check input files
[ -f "${INPUT_VCF}" ] || { echo "Error: Input VCF ${INPUT_VCF} not found."; exit 1; }
[ -f "${PED_FILE}" ] || { echo "Error: Pedigree file ${PED_FILE} not found."; exit 1; }
[ -f "${REFERENCE}" ] || { echo "Error: Reference ${REFERENCE} not found."; exit 1; }
[ -f "${REFERENCE}.fai" ] || { echo "Error: Reference index ${REFERENCE}.fai not found."; exit 1; }

# Run the annotation directly
echo "Annotating ${FAMILY} with PossibleDeNovo..."
gatk --java-options "-Xmx16g" VariantAnnotator \
  -R "${REFERENCE}" \
  -V "${INPUT_VCF}" \
  -O "${OUTPUT_VCF}" \
  -A PossibleDeNovo \
  --pedigree "${PED_FILE}"
[ $? -eq 0 ] || { echo "Error: Annotation failed for ${FAMILY}"; exit 1; }
echo "Annotation complete for ${FAMILY}"
