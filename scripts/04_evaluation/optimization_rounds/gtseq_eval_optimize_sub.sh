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
#SBATCH --nodelist=cbsusard

# --- Configuration ---

# Define project paths for robustness
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_2/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_3/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/349loci/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/399loci/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/299loci/eval/"
#PROJECT_ROOT="/local/workdir/cauret/gtseq_cisco/seq_panel_5/eval/"
#PROJECT_ROOT="/home/SardLab/projects/24/cisco_larvae_pedigree/seq_data/temp_folder/cisco_2024_test/R1/eval/"
PROJECT_ROOT="/home/SardLab/projects/26/chaumont_larvae_AA_gtseq_2025/TOW_7/eval/"
RSCRIPT_PATH="/home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_optimize.R"
#PRIMER_TAB="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/tm57/cisco_100wind_relax_primers_ordering_format_all_dimers.csv"
PRIMER_TAB="/home/SardLab/projects/24/cisco_larvae_pedigree/seq_data/temp_folder/cisco_2024_test/R1/AA_cisco_GTprimers.csv"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_2"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_3"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/349loci"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/399loci"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/299loci"
#FASTQ_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_5"
#FASTQ_DIR="/home/SardLab/projects/24/cisco_larvae_pedigree/seq_data/temp_folder/cisco_2024_test/R1"
FASTQ_DIR="/home/SardLab/projects/26/chaumont_larvae_AA_gtseq_2025/TOW_7"
OUT_DIR=${PROJECT_ROOT}
#TO_REMOVE_PRIMER_FILE="${PROJECT_ROOT}primers_to_remove_for_run2.csv"
#TO_REMOVE_PRIMER_FILE="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_to_remove_for_run2.csv"
#TO_REMOVE_PRIMER_FILE="/local/workdir/cauret/gtseq_cisco/seq_panel_3/primers_to_remove_for_runs2_3_combined.csv"
#TO_REMOVE_PRIMER_FILE="/local/workdir/cauret/gtseq_cisco/seq_panel_4/primers_to_remove_combined_loci_349.csv"
#TO_REMOVE_PRIMER_FILE="/local/workdir/cauret/gtseq_cisco/seq_panel_4/primers_to_remove_combined_loci_399.csv"
#TO_REMOVE_PRIMER_FILE="/local/workdir/cauret/gtseq_cisco/seq_panel_4/primers_to_remove_combined_loci_299.csv"
#TO_REMOVE_PRIMER_FILE="/local/workdir/cauret/gtseq_cisco/seq_panel_4/primers_to_remove_combined_loci_367.csv"
TO_REMOVE_PRIMER_FILE="/home/SardLab/projects/24/cisco_larvae_pedigree/seq_data/temp_folder/cisco_2024_test/R1/empty_to_remove.csv"
PRIMER_MISMATCH=3

#setting working directory
cd "$PROJECT_ROOT"

# --- Load R version ---
module load R/4.4.2

# --- Execution ---
# Access the array task ID within the script
echo "Processing task ID: $SLURM_ARRAY_TASK_ID"
#/usr/bin/time -v Rscript "$RSCRIPT_PATH" ${SLURM_ARRAY_TASK_ID}
/usr/bin/time -v Rscript "$RSCRIPT_PATH" ${SLURM_ARRAY_TASK_ID} ${PRIMER_TAB} ${FASTQ_DIR} ${OUT_DIR} ${TO_REMOVE_PRIMER_FILE} ${PRIMER_MISMATCH}

#sbatch --array=1-192%50 /home2/SardLab/cauret/scripts/gtseq_eval/gtseq_eval_optimize_sub.sh
