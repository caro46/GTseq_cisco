#!/bin/bash

#SBATCH --job-name=gtseq_eval_plot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=1G
#SBATCH --output=gtseq_eval_plot.%J.out
#SBATCH --error=gtseq_eval_plot.%J.err
#SBATCH --mail-user=caroline.tupper@oswego.edu
#SBATCH --mail-type=ALL

# --- Configuration ---

# Define project paths for robustness
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_2/eval/"
PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_3/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_2/eval/test/"
#RSCRIPT_PATH="/home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_optimize_plot.R"
RSCRIPT_PATH="/home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_optimize_plot_no_dil.R"
#OUT_DIR="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_2_files/"
OUT_DIR="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_3_files/"
F_CUMULATIVE_LIMIT="0.2"
#RUN_NUMBER=3
RUN_NUMBER=4

#setting working directory
cd "$PROJECT_ROOT"

# --- Load R version ---
module load R/4.4.2

# --- Execution ---
/usr/bin/time -v Rscript "$RSCRIPT_PATH" ${PROJECT_ROOT} ${OUT_DIR} ${F_CUMULATIVE_LIMIT} ${RUN_NUMBER}

#sbatch /home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_optimize_plot_sub.sh
