#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake-himem
#SBATCH --mem=12000

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"

[ -z "$FAMILY" ] || [ -z "$SAMPLE" ] && { echo "ERROR: FAMILY and SAMPLE must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[ -f "${BASE_DIR}/rg/${SAMPLE}.merged.bam" ] || { echo "ERROR: Input BAM not found: ${BASE_DIR}/rg/${SAMPLE}.merged.bam"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0

command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx10G" MarkDuplicates \
    --INPUT ${BASE_DIR}/rg/${SAMPLE}.merged.bam \
    --OUTPUT ${BASE_DIR}/${SAMPLE}.aligned.unsorted.dedup.bam \
    --METRICS_FILE ${BASE_DIR}/metrics/${SAMPLE}.duplicate_metrics \
    --VALIDATION_STRINGENCY SILENT \
    --ASSUME_SORT_ORDER "queryname" \
    --CLEAR_DT "false" \
    --READ_NAME_REGEX null \
    --CREATE_INDEX true

[ -f "${BASE_DIR}/${SAMPLE}.aligned.unsorted.dedup.bam" ] || { echo "ERROR: Output BAM not created: ${BASE_DIR}/${SAMPLE}.aligned.unsorted.dedup.bam"; exit 1; }
