#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH --mem=8000
#SBATCH -o logs/split_vcf_per_sample_%j.out
#SBATCH -e logs/split_vcf_per_sample_%j.err

FAMILY=$1
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"

[ -z "$FAMILY" ] && { echo "ERROR: FAMILY must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz" ] || { echo "ERROR: Filtered SNPs VCF not found"; exit 1; }
[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz" ] || { echo "ERROR: Filtered INDELs VCF not found"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load bcftools/1.19

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }

# Split SNPs
bcftools query -l ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz | while read SAMPLE; do
    bcftools view -s "$SAMPLE" ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz | bcftools view -i 'GT="alt"' -Oz -o ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz
    bcftools index ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.snps.vcf.gz
done

# Split INDELs
bcftools query -l ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz | while read SAMPLE; do
    bcftools view -s "$SAMPLE" ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz | bcftools view -i 'GT="alt"' -Oz -o ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz
    bcftools index ${BASE_DIR}/vcf/${SAMPLE}.hardfiltered.indels.vcf.gz
done
