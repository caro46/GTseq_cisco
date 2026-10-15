#!/bin/bash
#
# Example submission (default pattern: "_R1"):
#
# N=$(find /local/workdir/cauret/gtseq_cisco/seq_panel_5 \
#       -maxdepth 1 -type f \
#       \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" \) \
#     | wc -l)
# echo "$N"
# sbatch \
#   --array=1-${N} \
#   /home2/SardLab/cauret/scripts/bwa_amplicon.sh \
#   /local/workdir/cauret/gtseq_cisco/seq_panel_5 \
#   /local/workdir/cauret/gtseq_cisco/amplicon_info/cisco_amplicon.fa \
#   /local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons
#
#
# Example submission (custom pattern: ".R1."):
#
# N=$(find fastqs \
#       -maxdepth 1 -type f \
#       \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
#     | grep -F '.R1.' \
#     | wc -l)
#
# sbatch \
#   --array=1-${N} \
#   bwa_amplicon.sh \
#   fastqs \
#   amplicons.fa \
#   results \
#   .R1. \
#   .R2.
#

#SBATCH --job-name=bwa_amplicon
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=caroline.tupper@oswego.edu
#SBATCH --mail-type=ALL
#SBATCH --nodelist=cbsusard

set -euo pipefail

export PATH=/programs/bwa-mem2-2.3:$PATH

#########################################
# Inputs
#########################################

FASTQ_DIR=$1
REF=$2
OUTDIR=$3

# Optional, defaults to "_R1"
R1_PATTERN=${4:-_R1}

# Optional. Defaults to R1_PATTERN with "R1" swapped for "R2".
# If R1_PATTERN doesn't contain "R1", you MUST pass this explicitly.
R2_PATTERN=${5:-}
if [ -z "${R2_PATTERN}" ]; then
	if [[ "${R1_PATTERN}" != *R1* ]]; then
		echo "ERROR: R1_PATTERN ('${R1_PATTERN}') does not contain 'R1'," \
			"so R2_PATTERN cannot be inferred. Pass it explicitly as the 5th argument."
		exit 1
	fi
	R2_PATTERN=${R1_PATTERN/R1/R2}
fi

mkdir -p "${OUTDIR}"

#########################################
# Build list of R1 FASTQs (once, safely, shared across the whole array)
#########################################
#
# All array tasks race to get here at ~the same time. Only build the list
# once: task 1 writes it (via a temp file + atomic mv so partial reads are
# impossible), everyone else waits for the final file to show up.

FINAL_LIST="${OUTDIR}/r1_files.list"
TMP_LIST="${FINAL_LIST}.${SLURM_ARRAY_JOB_ID}.tmp"

if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
	find "${FASTQ_DIR}" \
		-maxdepth 1 \
		-type f \
		\( \
			-name "*${R1_PATTERN}*.fastq.gz" \
			-o \
			-name "*${R1_PATTERN}*.fq.gz" \
		\) \
	| sort \
	> "${TMP_LIST}"
	mv -f "${TMP_LIST}" "${FINAL_LIST}"
else
	WAITED=0
	while [ ! -f "${FINAL_LIST}" ]; do
		sleep 2
		WAITED=$((WAITED + 2))
		if [ "${WAITED}" -ge 300 ]; then
			echo "ERROR: Timed out waiting for ${FINAL_LIST} to be created by task 1"
			exit 1
		fi
	done
fi

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FINAL_LIST}")

if [ -z "${R1}" ]; then
	echo "No R1 file found for task ${SLURM_ARRAY_TASK_ID}"
	exit 1
fi

R2=${R1/${R1_PATTERN}/${R2_PATTERN}}

if [ ! -f "${R2}" ]; then
	echo "ERROR: Cannot find matching R2"
	echo "R1: ${R1}"
	echo "Expected R2: ${R2}"
	exit 1
fi

#########################################
# Sample name
#########################################

SAMPLE=$(basename "${R1}")
SAMPLE=${SAMPLE%%"${R1_PATTERN}"*}

echo "Processing ${SAMPLE}"

#########################################
# Verify bwa-mem2 index exists
#########################################

if [ ! -f "${REF}.bwt.2bit.64" ]; then
	echo "ERROR: bwa-mem2 index not found for:"
	echo "${REF}"
	echo
	echo "Build it once before submitting the array:"
	echo "bwa-mem2 index ${REF}"
	exit 1
fi

#########################################
# Read groups for GATK
#########################################

RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}"

#########################################
# Alignment (coordinate-sorted, GATK-ready)
#########################################

echo "Running bwa mem..."

bwa-mem2 mem \
	-M \
	-R "${RG}" \
	-t "${SLURM_CPUS_PER_TASK}" \
	"${REF}" \
	"${R1}" \
	"${R2}" \
| samtools sort \
	-@ "${SLURM_CPUS_PER_TASK}" \
	-m 1G \
	-o "${OUTDIR}/${SAMPLE}.bam"

# -M for GATK compatibility
samtools index "${OUTDIR}/${SAMPLE}.bam"

#########################################
# Keep primary proper pairs only
#########################################

samtools view \
	-b \
	-f 2 \
	-F 2304 \
	"${OUTDIR}/${SAMPLE}.bam" \
	> "${OUTDIR}/${SAMPLE}.proper_pairs.bam"

#########################################
# Reorder by read name for FASTQ export
#########################################

samtools collate \
	-o "${OUTDIR}/${SAMPLE}.proper_pairs.collated.bam" \
	"${OUTDIR}/${SAMPLE}.proper_pairs.bam"

#########################################
# Export mapped paired FASTQs
#########################################

samtools fastq \
	-1 "${OUTDIR}/${SAMPLE}.mapped_R1.fastq.gz" \
	-2 "${OUTDIR}/${SAMPLE}.mapped_R2.fastq.gz" \
	-0 /dev/null \
	-s /dev/null \
	-n \
	"${OUTDIR}/${SAMPLE}.proper_pairs.collated.bam"

#########################################
# Mapping statistics
#########################################

samtools flagstat \
	"${OUTDIR}/${SAMPLE}.bam" \
	> "${OUTDIR}/${SAMPLE}.flagstat.txt"

samtools idxstats \
	"${OUTDIR}/${SAMPLE}.bam" \
	> "${OUTDIR}/${SAMPLE}.idxstats.txt"

#########################################
# Validate paired FASTQs
#########################################

R1_COUNT=$(( $(zcat "${OUTDIR}/${SAMPLE}.mapped_R1.fastq.gz" | wc -l) / 4 ))
R2_COUNT=$(( $(zcat "${OUTDIR}/${SAMPLE}.mapped_R2.fastq.gz" | wc -l) / 4 ))

if [ "${R1_COUNT}" -ne "${R2_COUNT}" ]; then
	echo "WARNING: R1 (${R1_COUNT}) != R2 (${R2_COUNT})"
else
	echo "Paired FASTQs validated (${R1_COUNT} read pairs)"
fi

#########################################
# Cleanup intermediates
#########################################

rm -f "${OUTDIR}/${SAMPLE}.proper_pairs.bam"
rm -f "${OUTDIR}/${SAMPLE}.proper_pairs.collated.bam"

echo "Done."