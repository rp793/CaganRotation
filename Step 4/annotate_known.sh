#!/bin/bash
# Script: annotate_known.sh
# Purpose: Annotate VCFs with known canine variants and create index
# Usage: sbatch annotate_known.sh <FAMILY> <SAMPLE> (SAMPLE can be "trio" or sample ID)

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
#SBATCH -o ${BASE_DIR}/logs/annotate_known_%j.out
#SBATCH -e ${BASE_DIR}/logs/annotate_known_%j.err

if [ $# -ne 2 ]; then
  echo "Error: Please provide FAMILY and SAMPLE as arguments."
  echo "Usage: $0 <FAMILY> <SAMPLE>"
  exit 1
fi

KNOWN_VCF="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/variants/canfam3renamed.vcf.gz"

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load bcftools/1.19/gcc/3thngzmz

[ -f "$KNOWN_VCF" ] || { echo "Error: Known variants VCF not found: $KNOWN_VCF"; exit 1; }
[ -f "${KNOWN_VCF}.tbi" ] || { echo "Error: Index for ${KNOWN_VCF} not found."; exit 1; }

if [ "$SAMPLE" = "trio" ]; then
  INPUT_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.vcf.gz"
  ANNOTATED_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.known.vcf.gz"
else
  INPUT_VCF="${BASE_DIR}/vcf/annotated/${SAMPLE}.hardfiltered.snps.annotated.vcf.gz"
  ANNOTATED_VCF="${BASE_DIR}/vcf/annotated/${SAMPLE}.hardfiltered.snps.annotated.known.vcf.gz"
fi

[ -f "$INPUT_VCF" ] || { echo "Error: Input VCF not found: $INPUT_VCF"; exit 1; }

echo "Annotating $SAMPLE with known variant IDs from Ensembl..."
bcftools annotate \
  -a "$KNOWN_VCF" \
  -c ID \
  -Oz -o "$ANNOTATED_VCF" \
  "$INPUT_VCF"

[ $? -eq 0 ] || { echo "Error: Annotation failed for $SAMPLE"; exit 1; }

echo "Indexing annotated VCF for $SAMPLE..."
tabix -p vcf "$ANNOTATED_VCF"

[ $? -eq 0 ] || { echo "Error: Indexing failed for $ANNOTATED_VCF"; exit 1; }

echo "Annotation and indexing complete for $SAMPLE"
