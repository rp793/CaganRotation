#!/usr/bin/env bash
#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 03:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake  
#SBATCH -o logs/qc-Lane_%A.out

# Load modules
. /etc/profile.d/modules.sh                
module purge                               
module load rhel8/default-ccl              
module load ceuadmin/openjdk/17.0.8+7      
module load anaconda/3.2019-10             
module load ceuadmin/gatk/4.4.0.0 

# Check if FAMILY and SAMPLE are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: FAMILY and SAMPLE ID must be provided. Usage: $0 <family> <sample_id>" >&2
    exit 1
fi

FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
INPUT_BAM="${BASE_DIR}/rg/${SAMPLE}.merged.bam"  # Matches mergebam.sh output

# Check if input BAM file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file $INPUT_BAM does not exist" >&2
    exit 1
fi

# Create output directory if it doesn’t exist
mkdir -p "${BASE_DIR}/metrics"

# Run GATK command
gatk CollectMultipleMetrics \
	--INPUT "$INPUT_BAM" \
	--OUTPUT "${BASE_DIR}/metrics/${SAMPLE}.readgroup" \
	--ASSUME_SORTED true \
	--PROGRAM null \
	--PROGRAM CollectBaseDistributionByCycle \
	--PROGRAM CollectInsertSizeMetrics \
	--PROGRAM MeanQualityByCycle \
	--PROGRAM QualityScoreDistribution \
	--METRIC_ACCUMULATION_LEVEL null \
	--METRIC_ACCUMULATION_LEVEL ALL_READS
