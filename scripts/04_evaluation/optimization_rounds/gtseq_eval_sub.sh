#!/bin/bash

#SBATCH --job-name=gtseq_eval
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=1G
#SBATCH --output=gtseq_eval.%J.out
#SBATCH --error=gtseq_eval.%J.err
#SBATCH --mail-user=caroline.tupper@oswego.edu
#SBATCH --mail-type=ALL

# --- Configuration ---
# Define project paths for robustness
PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_1/eval/"
RSCRIPT_PATH="/home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval.R"
#setting working directory
cd "$PROJECT_ROOT"

# --- Load R version ---
module load R/4.4.2

# --- Execution ---
# Access the array task ID within the script
echo "Processing task ID: $SLURM_ARRAY_TASK_ID"
/usr/bin/time -v Rscript "$RSCRIPT_PATH" ${SLURM_ARRAY_TASK_ID}

#sbatch --array=1-188%50 /home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_sub.sh
