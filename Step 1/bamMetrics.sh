#!/usr/bin/env bash
#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time 12:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake 
#SBATCH -o logs/bamMetrics-%j.out
#SBATCH -e logs/bamMetrics-%j.err

# Load modules
. /etc/profile.d/modules.sh                
module purge                               
module load rhel8/default-ccl              
module load ceuadmin/openjdk/11.0.20+8  
module load anaconda/3.2019-10             
module load picard/2.9.2
module load samtools
export TMPDIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"

# Define Picard JAR path explicitly
PICARD_JAR="/usr/local/Cluster-Apps/picard/2.9.2/picard.jar"

# Check if Picard JAR exists
if [ ! -f "$PICARD_JAR" ]; then
    echo "Error: Picard JAR file $PICARD_JAR not found" >&2
    exit 1
fi
echo "Using Picard JAR: $PICARD_JAR"

# Check if FAMILY and SAMPLE are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: FAMILY and SAMPLE ID must be provided. Usage: $0 <family> <sample_id>" >&2
    exit 1
fi

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
METRICS_DIR="${BASE_DIR}/metrics"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"
BAM_FILE="${BASE_DIR}/${SAMPLE}-baserecal.bam"  # Updated to match applybqsr.sh
STATS_FILE="${METRICS_DIR}/${SAMPLE}.stats"

# Check if BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file $BAM_FILE does not exist" >&2
    exit 1
fi

# Check if reference file exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE does not exist" >&2
    exit 1
fi
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Error: Reference index ${REFERENCE}.fai does not exist" >&2
    exit 1
fi

# Create metrics directory
mkdir -p "$METRICS_DIR"
if [ ! -w "$METRICS_DIR" ]; then
    echo "Error: Cannot write to $METRICS_DIR" >&2
    exit 1
fi

# Run commands with error checking
echo "Running samtools flagstat..."
samtools flagstat "$BAM_FILE" > "${METRICS_DIR}/${SAMPLE}.flagstat" || {
    echo "Error: samtools flagstat failed" >&2
    exit 1
}

echo "Running picard CollectWgsMetrics..."
java -jar "$PICARD_JAR" CollectWgsMetrics INPUT="$BAM_FILE" O="${METRICS_DIR}/${SAMPLE}.wgs.txt" REFERENCE_SEQUENCE="$REFERENCE" STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT 2>> "${SLURM_JOB_ID}.picard_wgs.err"
if [ ! -f "${METRICS_DIR}/${SAMPLE}.wgs.txt" ]; then
    echo "Error: CollectWgsMetrics failed to produce output. Check ${SLURM_JOB_ID}.picard_wgs.err" >&2
    exit 1
fi

echo "Running picard CollectInsertSizeMetrics..."
java -jar "$PICARD_JAR" CollectInsertSizeMetrics INPUT="$BAM_FILE" O="${METRICS_DIR}/${SAMPLE}.insert_size.txt" H="${METRICS_DIR}/${SAMPLE}.insert_size_histogram.pdf" VALIDATION_STRINGENCY=LENIENT 2>> "${SLURM_JOB_ID}.picard_insert.err"
if [ ! -f "${METRICS_DIR}/${SAMPLE}.insert_size.txt" ]; then
    echo "Error: CollectInsertSizeMetrics failed to produce output. Check ${SLURM_JOB_ID}.picard_insert.err" >&2
    exit 1
fi

# Print to log
echo "$SAMPLE"
echo "..."
grep 'in total' "${METRICS_DIR}/${SAMPLE}.flagstat"
grep '%' "${METRICS_DIR}/${SAMPLE}.flagstat" | grep -v 'singleton'
echo "..."
head -8 "${METRICS_DIR}/${SAMPLE}.insert_size.txt" | tail -2 | cut -f 6,8
echo "..."
head -8 "${METRICS_DIR}/${SAMPLE}.wgs.txt" | tail -2 | cut -f 2,3,15,16,18,20-22
echo "..."

# Write to stats file
echo -e "BAM\t`stat -c '%y' "$BAM_FILE"`" >> "$STATS_FILE"
if [ -f "${BASE_DIR}/${SAMPLE}-baserecal.bam.md5" ]; then  # Updated to match applybqsr.sh
    md5sum --check "${BASE_DIR}/${SAMPLE}-baserecal.bam.md5" && echo -e "bamMd5\tVerified" >> "$STATS_FILE" || echo -e "bamMd5\tFailed" >> "$STATS_FILE"
fi
echo -e "totalReads\t`grep 'in total' "${METRICS_DIR}/${SAMPLE}.flagstat" | cut -f 1 -d ' '`" >> "$STATS_FILE"
echo -e "primaryReads\t`grep 'primary' "${METRICS_DIR}/${SAMPLE}.flagstat" | grep -v mapped | grep -v duplicates | cut -f 1 -d ' '`" >> "$STATS_FILE"
echo -e "mappedReads\t`grep 'mapped (' "${METRICS_DIR}/${SAMPLE}.flagstat" | grep -v 'primary' | cut -f 1 -d ' '`" >> "$STATS_FILE"
echo -e "primaryMapped\t`grep 'primary mapped' "${METRICS_DIR}/${SAMPLE}.flagstat" | cut -f 1 -d ' '`" >> "$STATS_FILE"
echo -e "pairedReads\t`grep 'properly paired' "${METRICS_DIR}/${SAMPLE}.flagstat" | cut -f 1 -d ' '`" >> "$STATS_FILE"
echo -e "meanInsertSize\t`head -8 "${METRICS_DIR}/${SAMPLE}.insert_size.txt" | tail -1 | cut -f 6`" >> "$STATS_FILE"
echo -e "meanReadDepth\t`head -8 "${METRICS_DIR}/${SAMPLE}.wgs.txt" | tail -1 | cut -f 2`" >> "$STATS_FILE"
