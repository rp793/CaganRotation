#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake-himem
#SBATCH --mem=12000
#SBATCH -o logs/filter_vcf_v2_%j.out
#SBATCH -e logs/filter_vcf_v2_%j.err

FAMILY=$1
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"

[ -z "$FAMILY" ] && { echo "ERROR: FAMILY must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.vcf.gz" ] || { echo "ERROR: Cohort VCF not found: ${BASE_DIR}/vcf/${FAMILY}_cohort.vcf.gz"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0

command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }

# Step 1: Split into SNPs and INDELs
gatk SelectVariants \
    -V ${BASE_DIR}/vcf/${FAMILY}_cohort.vcf.gz \
    -select-type SNP \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort.snps.vcf.gz

gatk SelectVariants \
    -V ${BASE_DIR}/vcf/${FAMILY}_cohort.vcf.gz \
    -select-type INDEL \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort.indels.vcf.gz

[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.snps.vcf.gz" ] || { echo "ERROR: SNPs VCF not created: ${BASE_DIR}/vcf/${FAMILY}_cohort.snps.vcf.gz"; exit 1; }
[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.indels.vcf.gz" ] || { echo "ERROR: INDELs VCF not created: ${BASE_DIR}/vcf/${FAMILY}_cohort.indels.vcf.gz"; exit 1; }

# Step 2: Hard filter SNPs
gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx10G" VariantFiltration \
    -R ${REFERENCE} \
    -V ${BASE_DIR}/vcf/${FAMILY}_cohort.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz

[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz" ] || { echo "ERROR: Hardfiltered SNPs VCF not created: ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.snps.vcf.gz"; exit 1; }

# Step 3: Hard filter INDELs
gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx10G" VariantFiltration \
    -R ${REFERENCE} \
    -V ${BASE_DIR}/vcf/${FAMILY}_cohort.indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz

[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz" ] || { echo "ERROR: Hardfiltered INDELs VCF not created: ${BASE_DIR}/vcf/${FAMILY}_cohort.hardfiltered.indels.vcf.gz"; exit 1; }
