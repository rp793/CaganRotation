#!/bin/bash
# Script: generate_dnm_input_v4.sh
# Purpose: Generate DNM input file for PhaseMyDeNovo using filtered de novo VCF
# Usage: sbatch generate_dnm_input.sh <FAMILY> <CHILD> <FATHER> <MOTHER>

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH -o /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}/logs/generate_dnm_%j.out
#SBATCH -e /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}/logs/generate_dnm_%j.err

FAMILY=$1
CHILD=$2
FATHER=$3
MOTHER=$4
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"

if [ $# -ne 4 ]; then
  echo "Error: Please provide FAMILY, CHILD, FATHER, and MOTHER as arguments."
  echo "Usage: $0 <FAMILY> <CHILD> <FATHER> <MOTHER>"
  exit 1
fi

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load bcftools/1.19/gcc/3thngzmz

DNM_INPUT="${BASE_DIR}/vcf/annotated/${CHILD}_dnm_input.txt"
FILTERED_VCF="${BASE_DIR}/vcf/annotated/${CHILD}.hardfiltered.snps.annotated.known.filtered.vcf.gz"
TRIO_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.vcf.gz"
BAM="${BASE_DIR}/${CHILD}-baserecal.bam"

# Check files
[ -f "$FILTERED_VCF" ] || { echo "Error: Filtered de novo VCF not found: $FILTERED_VCF"; exit 1; }
[ -f "$TRIO_VCF" ] || { echo "Error: Trio VCF not found: $TRIO_VCF"; exit 1; }
[ -f "$BAM" ] || { echo "Error: BAM not found: $BAM"; exit 1; }

echo "Generating DNM input for $CHILD (all calls)..."
echo -e "id\tchrom\tpos\tref\talt\tvcfs\tvcf_ids\tcram" > "$DNM_INPUT"
bcftools view -H "$FILTERED_VCF" | awk -F "\t" -v child="$CHILD" -v father="$FATHER" -v mother="$MOTHER" \
  '{print child"\t"$1"\t"$2"\t"$4"\t"$5"\t'"${TRIO_VCF}"'|'"${TRIO_VCF}"'|'"${TRIO_VCF}"'\t"child"|"father"|"mother"\t'"$BAM"'"}' >> "$DNM_INPUT"

[ $? -eq 0 ] && [ -s "$DNM_INPUT" ] || { echo "Error: Failed to generate DNM input for $CHILD"; exit 1; }

echo "DNM input generation complete for $CHILD"
