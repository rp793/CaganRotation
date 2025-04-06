#!/usr/bin/env bash
# Usage: sbatch -A ACCOUNT_NAME masterscript_vcf_v3.sh FAMILY_FOLDER
# Example: sbatch -A CAGAN-SL2-CPU masterscript_vcf_v3.sh GS-1

if [ $# -ne 1 ]; then
    echo "Usage: sbatch -A ACCOUNT_NAME $0 FAMILY_FOLDER"
    echo "Example: sbatch -A CAGAN-SL2-CPU $0 GS-1"
    exit 1
fi

FAMILY=$1
ACCOUNT=$SLURM_JOB_ACCOUNT
SCRIPT_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/scripts"
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"
REFERENCE="/home/rp793/rds/hpc-work/genomes/dog/Genomes/3.1/CanFam3.1_genomic.fna"

[ -z "$ACCOUNT" ] && { echo "ERROR: No account specified. Use -A ACCOUNT_NAME"; exit 1; }
[ -d "$BASE_DIR" ] || { echo "ERROR: Family directory not found: $BASE_DIR"; exit 1; }
[ -d "$SCRIPT_DIR" ] || { echo "ERROR: Script directory not found: $SCRIPT_DIR"; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: Reference genome not found: $REFERENCE"; exit 1; }

mkdir -p "${BASE_DIR}/logs" "${BASE_DIR}/gvcf" "${BASE_DIR}/vcf"
[ $? -eq 0 ] || { echo "ERROR: Failed to create subdirectories in $BASE_DIR"; exit 1; }

BAM_FILES=($(ls ${BASE_DIR}/*-baserecal.bam 2>/dev/null))
[ ${#BAM_FILES[@]} -eq 0 ] && { echo "ERROR: No *-baserecal.bam files found in $BASE_DIR"; exit 1; }

# Step 1: CombineGVCFs and GenotypeGVCFs per chromosome
jid_genotype=$(sbatch --parsable -A "${ACCOUNT}" -J "${FAMILY}.genotype_per_chr" "${SCRIPT_DIR}/genotype_gvcfs_per_chr.sh" "${FAMILY}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit genotype_gvcfs_per_chr job"; exit 1; }

# Step 2: Merge per-chromosome VCFs
jid_merge=$(sbatch --parsable -A "${ACCOUNT}" -J "${FAMILY}.merge_vcfs" --dependency=afterok:${jid_genotype} "${SCRIPT_DIR}/merge_vcfs.sh" "${FAMILY}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit merge_vcfs job"; exit 1; }

# Step 3: Filter and split into SNPs/INDELs
jid_filter=$(sbatch --parsable -A "${ACCOUNT}" -J "${FAMILY}.filter" --dependency=afterok:${jid_merge} "${SCRIPT_DIR}/filter_vcf_v3.sh" "${FAMILY}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit filter_vcf_v3 job"; exit 1; }

# Step 4: Split cohort VCFs into per-sample VCFs
jid_split=$(sbatch --parsable -A "${ACCOUNT}" -J "${FAMILY}.split" --dependency=afterok:${jid_filter} "${SCRIPT_DIR}/split_vcf_per_sample.sh" "${FAMILY}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit split_vcf_per_sample job"; exit 1; }

# Step 5: Generate per-sample VCF stats
# Submit without predefining array range; let vcf_stats_per_sample.sh handle it
jid_stats=$(sbatch --parsable -A "${ACCOUNT}" -J "${FAMILY}.stats" --dependency=afterok:${jid_split} "${SCRIPT_DIR}/vcf_stats_v3.sh" "${FAMILY}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit vcf_stats_per_sample job"; exit 1; }

echo "Submitted VCF pipeline jobs for family ${FAMILY}:"
echo "GenotypeGVCFsPerChr: ${jid_genotype}"
echo "MergeVCFs: ${jid_merge}"
echo "FilterVCF: ${jid_filter}"
echo "SplitPerSample: ${jid_split}"
echo "VCFStatsPerSample: ${jid_stats}"
