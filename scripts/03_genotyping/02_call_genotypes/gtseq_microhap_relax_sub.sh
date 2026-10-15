#!/bin/bash

#SBATCH --job-name=gtseq_microhap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=20G
#SBATCH --output=gtseq_microhap.%J.out
#SBATCH --error=gtseq_microhap.%J.err
#SBATCH --mail-user=caroline.tupper@oswego.edu
#SBATCH --mail-type=ALL
#SBATCH --nodelist=cbsusard

# --- Configuration ---

# Path to program
SIF_PATH="/home2/SardLab/cauret/software/gtseq_microhap"
BWA_PATH="/programs/bwa-0.7.17/bwa"
PRIMERS="$1"
FASTQ="$2" # Directory containing R1/R2 fastqs
OUTDIR="$3"
REF="$4"
MAX_AMPLICON="${5:-150}" # Use arg 5, default to 150
VCF_MIN_PID=0.8
MIN_DP=5

# --- Execution ---

apptainer run --bind /home --bind /home2 --bind /local --bind /programs "$SIF_PATH"/gtseq.sif --primers ${PRIMERS} \
	--indir ${FASTQ} \
	--outdir ${OUTDIR} \
	--reference-fasta ${REF} \
	--max-amplicon ${MAX_AMPLICON} --emit-vcf --vcf-build-bwa-index --bwa-path ${BWA_PATH} \
	--vcf-min-pid ${VCF_MIN_PID} --min-depth ${MIN_DP}

# --- Usage ---
# Uses default max amplicon of 150
#sbatch /home2/SardLab/cauret/scripts/gtseq_microhap_sub.sh /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_5_files/locus_primers_367.csv /local/workdir/cauret/gtseq_cisco/seq_panel_5/ /local/workdir/cauret/gtseq_cisco/seq_panel_5/gtseq_microhap /home2/SardLab/cauret/reference_genomes/Coregonus_artedi_v1.fa