
```bash
#!/usr/bin/env bash

# Check if both family folder and sample ID are provided
if [ $# -ne 2 ]; then
    echo "Usage: sbatch -A ACCOUNT_NAME masterscript.sh FAMILY_FOLDER SAMPLE_ID"
    echo "Example: sbatch -A CAGAN-SL2-CPU masterscript.sh GS-1 SAMPLE123"
    exit 1
fi

FAMILY=$1
SAMPLE=$2
ACCOUNT=$SLURM_JOB_ACCOUNT
SCRIPT_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/scripts"
BASE_DIR="/home/rp793/rds/rds-dog-de-novo-mXW0tnsCK5M/data/dogtest/${FAMILY}"

# Check if ACCOUNT is set
if [ -z "$ACCOUNT" ]; then
    echo "ERROR: No account specified. Use -A ACCOUNT_NAME with sbatch."
    exit 1
fi

# Check if BASE_DIR exists
if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Family directory does not exist: $BASE_DIR"
    exit 1
fi

# Create family-specific subdirectories for outputs
mkdir -p "${BASE_DIR}/logs" "${BASE_DIR}/metrics" "${BASE_DIR}/rg"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create subdirectories in $BASE_DIR"
    exit 1
fi

# Check if SCRIPT_DIR exists and is readable
if [ ! -d "$SCRIPT_DIR" ] || [ ! -r "$SCRIPT_DIR" ]; then
    echo "ERROR: Script directory not found or not readable: $SCRIPT_DIR"
    exit 1
fi

# Submit main pipeline jobs
jid1=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.fastq2sam" "${SCRIPT_DIR}/fastq2sam.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit fastq2sam job"; exit 1; }
jid2=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.markadapters" --dependency=afterok:${jid1} "${SCRIPT_DIR}/markadapters.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit markadapters job"; exit 1; }
jid3=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.sam2fastq" --dependency=afterok:${jid2} "${SCRIPT_DIR}/sam2fastq.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit sam2fastq job"; exit 1; }
jid4=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.alignfastq" --dependency=afterok:${jid3} "${SCRIPT_DIR}/alignfastq.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit alignfastq job"; exit 1; }
jid5=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.mergebam" --dependency=afterok:${jid4} "${SCRIPT_DIR}/mergebam.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit mergebam job"; exit 1; }
jid6=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.markdupl" --dependency=afterok:${jid5} "${SCRIPT_DIR}/markdupl.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit markdupl job"; exit 1; }
jid7=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.sortbam" --dependency=afterok:${jid6} "${SCRIPT_DIR}/sortbam.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit sortbam job"; exit 1; }
jid8=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.baserecal" --dependency=afterok:${jid7} "${SCRIPT_DIR}/baserecal.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit baserecal job"; exit 1; }
jid9=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.applybqsr" --dependency=afterok:${jid8} "${SCRIPT_DIR}/applybqsr.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit applybqsr job"; exit 1; }

# Submit metrics jobs after applybqsr completes
jid10=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.qcyield" --dependency=afterok:${jid1} "${SCRIPT_DIR}/qcyieldmetrics.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit qcyieldmetrics job"; exit 1; }
jid11=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.qclane" --dependency=afterok:${jid5} "${SCRIPT_DIR}/qclanergmetrics.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit qclanergmetrics job"; exit 1; }
jid12=$(sbatch --parsable -A "${ACCOUNT}" -J "${SAMPLE}.bammetrics" --dependency=afterok:${jid9} "${SCRIPT_DIR}/bamMetrics.sh" "${FAMILY}" "${SAMPLE}")
[ $? -eq 0 ] || { echo "ERROR: Failed to submit bamMetrics job"; exit 1; }

echo "Submitted jobs for sample ${SAMPLE} from family ${FAMILY}:"
echo "Fastq2Sam: ${jid1}"
echo "MarkAdapters: ${jid2}"
echo "Sam2Fastq: ${jid3}"
echo "AlignFastq: ${jid4}"
echo "MergeBam: ${jid5}"
echo "MarkDuplicates: ${jid6}"
echo "SortBam: ${jid7}"
echo "BaseRecalibrator: ${jid8}"
echo "ApplyBQSR: ${jid9}"
echo "QualityYieldMetrics: ${jid10}"
echo "LaneMetrics: ${jid11}"
echo "BamMetrics: ${jid12}"
```
