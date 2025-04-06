#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH -o logs/vcf_stats_per_sample_%A_%a.out
#SBATCH -e logs/vcf_stats_per_sample_%A_%a.err
# Removed #SBATCH --array=0-%N since master script no longer sets it

FAMILY=$1
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"

[ -z "$FAMILY" ] && { echo "ERROR: FAMILY must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }

# Load modules
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load bcftools/1.19
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }

# Get list of samples from cohort VCF
COHORT_VCF="${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz"
[ -f "$COHORT_VCF" ] || { echo "ERROR: Cohort VCF not found: $COHORT_VCF"; exit 1; }
SAMPLES=($(bcftools query -l "$COHORT_VCF"))
[ ${#SAMPLES[@]} -eq 0 ] && { echo "ERROR: No samples found in cohort VCF"; exit 1; }

# If array job, use task ID; otherwise, process all samples
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
    [ -z "$SAMPLE" ] && { echo "ERROR: No sample for array task ID $SLURM_ARRAY_TASK_ID"; exit 1; }
else
    # Non-array mode: process all samples
    for SAMPLE in "${SAMPLES[@]}"; do
        STATS_FILE="${BASE_DIR}/vcf/${SAMPLE}.stats"
        [ -f "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz" ] || { echo "ERROR: Filtered SNPs VCF not found: ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz"; exit 1; }
        [ -f "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz" ] || { echo "ERROR: Filtered INDELs VCF not found: ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz"; exit 1; }

        bcftools stats -s- -f PASS "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz" > "${BASE_DIR}/vcf/${SAMPLE}_snps.vchk"
        bcftools stats -s- -f PASS "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz" > "${BASE_DIR}/vcf/${SAMPLE}_indels.vchk"

        echo -e "SNPs\t$(grep 'number of SNPs:' ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | cut -f 4)" > "$STATS_FILE"
        echo -e "snpHom\t$(grep PSC ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | tail -1 | cut -f 5)" >> "$STATS_FILE"
        echo -e "snpHet\t$(grep PSC ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | tail -1 | cut -f 6)" >> "$STATS_FILE"
        echo -e "InDels\t$(grep 'number of indels:' ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | cut -f 4)" >> "$STATS_FILE"
        echo -e "insHet\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 8)" >> "$STATS_FILE"
        echo -e "delHet\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 9)" >> "$STATS_FILE"
        echo -e "insHom\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 10)" >> "$STATS_FILE"
        echo -e "delHom\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 11)" >> "$STATS_FILE"
        echo -e "dateCompleted\t$(date +"%Y-%m-%d %H:%M:%S")" >> "$STATS_FILE"

        [ -s "$STATS_FILE" ] || { echo "ERROR: Stats file is empty: $STATS_FILE"; exit 1; }
    done
    exit 0
fi

# Array mode: process one sample
STATS_FILE="${BASE_DIR}/vcf/${SAMPLE}.stats"
[ -f "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz" ] || { echo "ERROR: Filtered SNPs VCF not found: ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz"; exit 1; }
[ -f "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz" ] || { echo "ERROR: Filtered INDELs VCF not found: ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz"; exit 1; }

bcftools stats -s- -f PASS "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz" > "${BASE_DIR}/vcf/${SAMPLE}_snps.vchk"
bcftools stats -s- -f PASS "${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz" > "${BASE_DIR}/vcf/${SAMPLE}_indels.vchk"

echo -e "SNPs\t$(grep 'number of SNPs:' ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | cut -f 4)" > "$STATS_FILE"
echo -e "snpHom\t$(grep PSC ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | tail -1 | cut -f 5)" >> "$STATS_FILE"
echo -e "snpHet\t$(grep PSC ${BASE_DIR}/vcf/${SAMPLE}_snps.vchk | tail -1 | cut -f 6)" >> "$STATS_FILE"
echo -e "InDels\t$(grep 'number of indels:' ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | cut -f 4)" >> "$STATS_FILE"
echo -e "insHet\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 8)" >> "$STATS_FILE"
echo -e "delHet\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 9)" >> "$STATS_FILE"
echo -e "insHom\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 10)" >> "$STATS_FILE"
echo -e "delHom\t$(grep PSI ${BASE_DIR}/vcf/${SAMPLE}_indels.vchk | tail -1 | cut -f 11)" >> "$STATS_FILE"
echo -e "dateCompleted\t$(date +"%Y-%m-%d %H:%M:%S")" >> "$STATS_FILE"

[ -s "$STATS_FILE" ] || { echo "ERROR: Stats file is empty: $STATS_FILE"; exit 1; }
