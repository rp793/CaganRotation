#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH --mem=16000
#SBATCH --array=0-38
#SBATCH -o logs/genotype_gvcfs_per_chr_%A_%a.out
#SBATCH -e logs/genotype_gvcfs_per_chr_%A_%a.err

FAMILY=$1
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"
CHROMOSOMES=(
    "NC_006583.3" "NC_006584.3" "NC_006585.3" "NC_006586.3" "NC_006587.3"
    "NC_006588.3" "NC_006589.3" "NC_006590.3" "NC_006591.3" "NC_006592.3"
    "NC_006593.3" "NC_006594.3" "NC_006595.3" "NC_006596.3" "NC_006597.3"
    "NC_006598.3" "NC_006599.3" "NC_006600.3" "NC_006601.3" "NC_006602.3"
    "NC_006603.3" "NC_006604.3" "NC_006605.3" "NC_006606.3" "NC_006607.3"
    "NC_006608.3" "NC_006609.3" "NC_006610.3" "NC_006611.3" "NC_006612.3"
    "NC_006613.3" "NC_006614.3" "NC_006615.3" "NC_006616.3" "NC_006617.3"
    "NC_006618.3" "NC_006619.3" "NC_006620.3" "NC_006621.3"
)

CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID]}
[ -z "$FAMILY" ] && { echo "ERROR: FAMILY must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

# Gather per-chromosome gVCFs from haplotypeCallersplit.sh
GVCF_ARGS=""
for SAMPLE in $(ls ${BASE_DIR}/*-baserecal.bam | sed 's/.*\///; s/-baserecal.bam//'); do
    GVCF="${BASE_DIR}/gvcf/${SAMPLE}-chr_${CHR}.g.vcf.gz"  # Updated path
    [ -f "$GVCF" ] || { echo "ERROR: gVCF not found: $GVCF"; exit 1; }
    GVCF_ARGS+=" -V $GVCF"
done

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0

command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx14G" CombineGVCFs \
    -R ${REFERENCE} \
    ${GVCF_ARGS} \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort_chr_${CHR}.g.vcf.gz

[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort_chr_${CHR}.g.vcf.gz" ] || { echo "ERROR: Cohort gVCF not created for chr $CHR"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx14G" GenotypeGVCFs \
    -R ${REFERENCE} \
    -V ${BASE_DIR}/vcf/${FAMILY}_cohort_chr_${CHR}.g.vcf.gz \
    -O ${BASE_DIR}/vcf/${FAMILY}_cohort_chr_${CHR}.vcf.gz

[ -f "${BASE_DIR}/vcf/${FAMILY}_cohort_chr_${CHR}.vcf.gz" ] || { echo "ERROR: Cohort VCF not created for chr $CHR"; exit 1; }
