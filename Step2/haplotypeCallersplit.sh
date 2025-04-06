#!/usr/bin/env bash

#! RUN : sbatch haplotypeCallersplit.sh <FAMILY> <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # For multi-threading
#SBATCH --mem=16G          # Matches -Xmx16G
#SBATCH --time 12:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=rp793@cam.ac.uk
#SBATCH -p cclake
#SBATCH --account=CAGAN-SL2-CPU 
#SBATCH --array=0-38
#SBATCH -o logs/haplotypeCallersplit-%A_%a.out

# Check if both FAMILY and SAMPLE are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: FAMILY and SAMPLE ID must be provided." >&2
    echo "Usage: sbatch $0 <family> <sample_id>" >&2
    echo "Example: sbatch $0 GS-1 SAMPLE123" >&2
    exit 1
fi

# Load modules
. /etc/profile.d/modules.sh                
module purge                               
module load rhel8/default-ccl              
module load ceuadmin/openjdk/17.0.8+7      
module load anaconda/3.2019-10             
module load ceuadmin/gatk/4.4.0.0 

# Set variables
FAMILY=$1
SAMPLE=$2
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
TMP_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/tmp"
METRICS_DIR="${BASE_DIR}/metrics"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"
PCR_MODEL='CONSERVATIVE'

# Check critical files
if [ ! -f "$REFERENCE" ] || [ ! -f "${REFERENCE}.fai" ] || [ ! -f "${REFERENCE%.fna}.dict" ]; then
    echo "Reference file or indices missing" >&2
    exit 1
fi
if [ ! -f "${BASE_DIR}/${SAMPLE}-baserecal.bam" ]; then
    echo "BAM file missing: ${BASE_DIR}/${SAMPLE}-baserecal.bam" >&2
    exit 1
fi

# Array of chromosome names
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

# Select chromosome based on SLURM array task ID
CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID]}
OUTPUT="${BASE_DIR}/gvcf/${SAMPLE}-chr_${CHR}.g.vcf.gz"  # Updated to include BASE_DIR and gvcf subdirectory

# Ensure gvcf directory exists
mkdir -p "${BASE_DIR}/gvcf"

echo "Processing sample ${SAMPLE} for chromosome ${CHR} in family ${FAMILY}"

gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} -Xmx16G" HaplotypeCaller \
  -R ${REFERENCE} \
  -I ${BASE_DIR}/${SAMPLE}-baserecal.bam \
  -O ${OUTPUT} \
  -L ${CHR} \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  -ERC GVCF \
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
  --pcr-indel-model ${PCR_MODEL} \
  --native-pair-hmm-threads 4

if [ $? -ne 0 ]; then
    echo "HaplotypeCaller failed for sample ${SAMPLE} on chromosome ${CHR} in family ${FAMILY}" >> logs/haplotypeCallersplit-${SAMPLE}.error.log
    exit 1
fi
