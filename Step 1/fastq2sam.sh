#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake-himem

set -e
FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"

# Check arguments
[ -z "$FAMILY" ] || [ -z "$SAMPLE" ] && { echo "ERROR: FAMILY and SAMPLE must be provided"; exit 1; }

# Check directories
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0

# Check if gatk is available
command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }

FILES=(
    "${BASE_DIR}/${SAMPLE}_1.fastq.gz"
    "${BASE_DIR}/${SAMPLE}_2.fastq.gz"
)

for file in "${FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Input file not found: $file"
        exit 1
    fi
done

FASTQ=$(zcat ${FILES[0]} | head -1)
[ -z "$FASTQ" ] && { echo "ERROR: Failed to read FASTQ header from ${FILES[0]}"; exit 1; }

LIBRARY=${SAMPLE}
FLOWCELL=$(echo "$FASTQ" | cut -f 3 -d':')
LANE=$(echo "$FASTQ" | cut -f 4 -d':')
[ -z "$FLOWCELL" ] || [ -z "$LANE" ] && { echo "ERROR: Failed to extract FLOWCELL or LANE from FASTQ header"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx4G" \
    FastqToSam \
    --FASTQ ${FILES[0]} \
    --FASTQ2 ${FILES[1]} \
    --OUTPUT ${BASE_DIR}/rg/${SAMPLE}.unaligned.bam \
    --READ_GROUP_NAME ${FLOWCELL}.L${LANE} \
    --SAMPLE_NAME ${SAMPLE} \
    --LIBRARY_NAME ${LIBRARY} \
    --PLATFORM_UNIT ${FLOWCELL}.L${LANE}.${SAMPLE} \
    --PLATFORM illumina \
    --SEQUENCING_CENTER CRUK-CI

[ -f "${BASE_DIR}/rg/${SAMPLE}.unaligned.bam" ] || { echo "ERROR: Output BAM not created: ${BASE_DIR}/rg/${SAMPLE}.unaligned.bam"; exit 1; }
