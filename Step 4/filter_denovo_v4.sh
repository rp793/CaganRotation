#!/bin/bash
# Script: filter_denovo_v4.sh
# Purpose: Annotate trio VCF with VAF and filter for de novo variants per child
# Usage: sbatch filter_denovo.sh <FAMILY> <CHILD>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH -o /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/%x/logs/filter_denovo_%j.out
#SBATCH -e /home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/%x/logs/filter_denovo_%j.err

FAMILY=$1
CHILD=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
PED_FILE="${BASE_DIR}/pedigree/${FAMILY}.ped"

if [ $# -ne 2 ]; then
  echo "Error: Please provide FAMILY and CHILD as arguments."
  echo "Usage: $0 <FAMILY> <CHILD>"
  exit 1
fi

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load bcftools/1.19/gcc/3thngzmz

CHROMOSOMES="NC_006583.3,NC_006584.3,NC_006585.3,NC_006586.3,NC_006587.3,NC_006588.3,NC_006589.3,NC_006590.3,NC_006591.3,NC_006592.3,NC_006593.3,NC_006594.3,NC_006595.3,NC_006596.3,NC_006597.3,NC_006598.3,NC_006599.3,NC_006600.3,NC_006601.3,NC_006602.3,NC_006603.3,NC_006604.3,NC_006605.3,NC_006606.3,NC_006607.3,NC_006608.3,NC_006609.3,NC_006610.3,NC_006611.3,NC_006612.3,NC_006613.3,NC_006614.3,NC_006615.3,NC_006616.3,NC_006617.3,NC_006618.3,NC_006619.3,NC_006620.3,NC_006621.3"
INPUT_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.known.vcf.gz"
TEMP_VCF="${BASE_DIR}/vcf/annotated/${FAMILY}.trio.hardfiltered.snps.annotated.known.vaf.vcf.gz"
FILTERED_VCF="${BASE_DIR}/vcf/annotated/${CHILD}.hardfiltered.snps.annotated.known.filtered.vcf.gz"

[ -f "$INPUT_VCF" ] || { echo "Error: Input VCF not found: $INPUT_VCF"; exit 1; }
[ -f "$PED_FILE" ] || { echo "Error: Pedigree file not found: $PED_FILE"; exit 1; }

# Determine sample indices
SAMPLES=($(bcftools query -l "$INPUT_VCF"))
CHILD_IDX=-1
FATHER_IDX=-1
MOTHER_IDX=-1
FATHER=$(awk '$3 == "0" && $5 == "1" {print $2}' "$PED_FILE" | head -n 1)
MOTHER=$(awk '$3 == "0" && $5 == "2" {print $2}' "$PED_FILE" | head -n 1)
for i in "${!SAMPLES[@]}"; do
  if [ "${SAMPLES[$i]}" = "$CHILD" ]; then CHILD_IDX=$i; fi
  if [ "${SAMPLES[$i]}" = "$FATHER" ]; then FATHER_IDX=$i; fi
  if [ "${SAMPLES[$i]}" = "$MOTHER" ]; then MOTHER_IDX=$i; fi
done
[ "$CHILD_IDX" -eq -1 ] || [ "$FATHER_IDX" -eq -1 ] || [ "$MOTHER_IDX" -eq -1 ] && { echo "Error: Could not map samples in VCF"; exit 1; }

echo "Annotating $FAMILY VCF with VAF..."
bcftools +fill-tags "$INPUT_VCF" -o "$TEMP_VCF.uncompressed" -- -t VAF 2> "${BASE_DIR}/logs/vaf_error_%j.log"
[ $? -eq 0 ] && [ -s "$TEMP_VCF.uncompressed" ] || { echo "Error: Failed to annotate VAF for $FAMILY, see vaf_error_%j.log"; exit 1; }
echo "VAF annotation completed, uncompressed file size: $(ls -lh $TEMP_VCF.uncompressed | awk '{print $5}')"

[ -f "$TEMP_VCF" ] && rm -f "$TEMP_VCF"
/usr/local/software/spack/csd3/opt-2024-06-01/linux-rocky8-cascadelake/gcc-13.3.0/htslib-1.19.1-ucjgd75np23amq72cviwr4iqolfzpgi6/bin/bgzip -c "$TEMP_VCF.uncompressed" > "$TEMP_VCF"
[ $? -eq 0 ] && [ -s "$TEMP_VCF" ] || { echo "Error: Compression failed for $TEMP_VCF"; ls -lh "$TEMP_VCF.uncompressed"; exit 1; }
echo "Compression completed, compressed file size: $(ls -lh $TEMP_VCF | awk '{print $5}')"

/usr/local/software/spack/csd3/opt-2024-06-01/linux-rocky8-cascadelake/gcc-13.3.0/htslib-1.19.1-ucjgd75np23amq72cviwr4iqolfzpgi6/bin/tabix -p vcf "$TEMP_VCF"
[ $? -eq 0 ] || { echo "Error: Failed to index $TEMP_VCF"; ls -lh "$TEMP_VCF"; exit 1; }
echo "Indexing completed"

rm -f "$TEMP_VCF.uncompressed"

echo "Filtering $CHILD for de novo variants..."
echo "CHILD=$CHILD (idx=$CHILD_IDX), FATHER=$FATHER (idx=$FATHER_IDX), MOTHER=$MOTHER (idx=$MOTHER_IDX)"
echo "VCF samples: ${SAMPLES[@]}"
FILTER_EXPR="FMT/VAF[${CHILD_IDX}:0] >= 0.3 && FMT/VAF[${CHILD_IDX}:0] <= 0.7 && FMT/VAF[${FATHER_IDX}:0] < 0.1 && FMT/VAF[${MOTHER_IDX}:0] < 0.1 && FMT/AD[${FATHER_IDX}:1] <= 2 && FMT/AD[${MOTHER_IDX}:1] <= 2 && FMT/DP[${CHILD_IDX}] > 8 && FMT/DP[${FATHER_IDX}] > 8 && FMT/DP[${MOTHER_IDX}] > 8 && INFO/hiConfDeNovo != \".\" && ID = \".\""
bcftools view \
  -i "$FILTER_EXPR" \
  -r "$CHROMOSOMES" \
  -s "$CHILD" \
  -Oz -o "$FILTERED_VCF" \
  "$TEMP_VCF"
  
[ $? -eq 0 ] || { echo "Error: Filtering failed for $CHILD"; exit 1; }
echo "Cleaning up temporary VCF..."
rm -f "$TEMP_VCF"
echo "Filtering complete for $CHILD"
