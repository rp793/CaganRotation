#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"

[ -z "$FAMILY" ] || [ -z "$SAMPLE" ] && { echo "ERROR: FAMILY and SAMPLE must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[ -f "${BASE_DIR}/${SAMPLE}.sorted.bam" ] || { echo "ERROR: Input BAM not found: ${BASE_DIR}/${SAMPLE}.sorted.bam"; exit 1; }
[ -f "${BASE_DIR}/${SAMPLE}.bsqr.txt" ] || { echo "ERROR: Recalibration table not found: ${BASE_DIR}/${SAMPLE}.bsqr.txt"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0

command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx16G" ApplyBQSR \
    -R ${REFERENCE} \
    -I ${BASE_DIR}/${SAMPLE}.sorted.bam \
    -O ${BASE_DIR}/${SAMPLE}-baserecal.bam \
    -bqsr ${BASE_DIR}/${SAMPLE}.bsqr.txt \
    --static-quantized-quals 10 \
    --static-quantized-quals 20 \
    --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities

[ -f "${BASE_DIR}/${SAMPLE}-baserecal.bam" ] || { echo "ERROR: Output BAM not created: ${BASE_DIR}/${SAMPLE}-baserecal.bam"; exit 1; }
