#!/usr/bin/env bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
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
[ -f "${BASE_DIR}/rg/${SAMPLE}.fastq.gz" ] || { echo "ERROR: Input FASTQ not found: ${BASE_DIR}/rg/${SAMPLE}.fastq.gz"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load ceuadmin/openjdk/17.0.8+7
module load anaconda/3.2019-10
module load ceuadmin/gatk/4.4.0.0
module load bwa/0.7.12
module load samtools

command -v bwa >/dev/null 2>&1 || { echo "ERROR: BWA not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: Samtools not found"; exit 1; }

bwa mem -K 100000000 -p -v 3 -t 16 -Y ${REFERENCE} ${BASE_DIR}/rg/${SAMPLE}.fastq.gz | samtools view -h -b -o ${BASE_DIR}/rg/${SAMPLE}.aligned.bam
[ $? -eq 0 ] || { echo "ERROR: BWA or samtools pipeline failed"; exit 1; }
[ -f "${BASE_DIR}/rg/${SAMPLE}.aligned.bam" ] || { echo "ERROR: Output BAM not created: ${BASE_DIR}/rg/${SAMPLE}.aligned.bam"; exit 1; }
