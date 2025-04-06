#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake-himem
#SBATCH --mem=25000

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"

[ -z "$FAMILY" ] || [ -z "$SAMPLE" ] && { echo "ERROR: FAMILY and SAMPLE must be provided"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Base directory not found: $BASE_DIR"; exit 1; }
[ -d "$TMP_DIR" ] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[ -f "${BASE_DIR}/rg/${SAMPLE}.unaligned.bam" ] || { echo "ERROR: Unmapped BAM not found: ${BASE_DIR}/rg/${SAMPLE}.unaligned.bam"; exit 1; }
[ -f "${BASE_DIR}/rg/${SAMPLE}.aligned.bam" ] || { echo "ERROR: Aligned BAM not found: ${BASE_DIR}/rg/${SAMPLE}.aligned.bam"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0
module load bwa/0.7.12
module load samtools

command -v gatk >/dev/null 2>&1 || { echo "ERROR: GATK not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: Samtools not found"; exit 1; }

PG_ID=$(samtools view -H ${BASE_DIR}/rg/${SAMPLE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 2 | sed 's/ID://')
PG_PN=$(samtools view -H ${BASE_DIR}/rg/${SAMPLE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 3 | sed 's/PN://')
PG_VN=$(samtools view -H ${BASE_DIR}/rg/${SAMPLE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 4 | sed 's/VN://')
PG_CL=$(samtools view -H ${BASE_DIR}/rg/${SAMPLE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 5 | sed 's/CL://')
[ -z "$PG_ID" ] && { echo "ERROR: Failed to extract PG_ID from BAM header"; exit 1; }

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx4G" MergeBamAlignment \
    --REFERENCE_SEQUENCE ${REFERENCE} \
    --UNMAPPED_BAM ${BASE_DIR}/rg/${SAMPLE}.unaligned.bam \
    --ALIGNED_BAM ${BASE_DIR}/rg/${SAMPLE}.aligned.bam \
    --OUTPUT ${BASE_DIR}/rg/${SAMPLE}.merged.bam \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --SORT_ORDER "unsorted" \
    --PROGRAM_RECORD_ID "${PG_ID}" \
    --PROGRAM_GROUP_VERSION "${PG_VN}" \
    --PROGRAM_GROUP_COMMAND_LINE "${PG_CL}" \
    --PROGRAM_GROUP_NAME "${PG_PN}" \
    --CLIP_ADAPTERS false \
    --UNMAP_CONTAMINANT_READS true \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --CREATE_INDEX false \
    --ADD_MATE_CIGAR true \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --INCLUDE_SECONDARY_ALIGNMENTS true

[ -f "${BASE_DIR}/rg/${SAMPLE}.merged.bam" ] || { echo "ERROR: Output BAM not created: ${BASE_DIR}/rg/${SAMPLE}.merged.bam"; exit 1; }
