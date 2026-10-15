#!/bin/bash
#SBATCH --job-name=rename_files_issues_GATK
#SBATCH --mail-user=caroline.tupper@oswego.edu
#SBATCH --mail-type=ALL
#SBATCH --nodelist=cbsusard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=rename_files_issues_GATK_%j.out
#SBATCH --error=rename_files_issues_GATK_%j.err

set -euo pipefail

##########################################################################
# EDIT THESE
##########################################################################

FASTA="/local/workdir/cauret/gtseq_cisco/amplicon_info/cisco_amplicon.fa"

BAM_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons"

OUT_DIR="/local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons_gatk"

##########################################################################

mkdir -p "$OUT_DIR"

OUT_FASTA="${OUT_DIR}/$(basename "${FASTA%.*}")_gatk.fa"

echo "===================================================="
echo "Creating GATK-safe FASTA"
echo "===================================================="

awk '
/^>/{
	sub(/\[.*\]/,"",$0)
}
{print}
' "$FASTA" > "$OUT_FASTA"

echo
echo "Checking for duplicate locus IDs..."

DUPS=$(grep '^>' "$OUT_FASTA" | sed 's/^>//' | sort | uniq -d)

if [[ -n "$DUPS" ]]; then
	echo "ERROR: duplicate contig names generated:"
	echo "$DUPS"
	exit 1
fi

echo "No duplicate contig IDs found."

echo
echo "===================================================="
echo "Reheadering BAMs"
echo "===================================================="

for BAM in "$BAM_DIR"/*.bam
do
	BASE=$(basename "$BAM" .bam)

	echo
	echo "Processing: $BASE"

	TMP_HEADER=$(mktemp)

	samtools view -H "$BAM" | \
	sed -E '/^@SQ/s/(SN:[^[]+)\[[^]]+\]/\1/' \
	> "$TMP_HEADER"

	OUT_BAM="${OUT_DIR}/${BASE}.bam"

	samtools reheader "$TMP_HEADER" "$BAM" > "$OUT_BAM"

	samtools index "$OUT_BAM"

	rm -f "$TMP_HEADER"

	echo "Finished: $BASE"
done

echo
echo "===================================================="
echo "DONE"
echo "===================================================="
echo "New FASTA:"
echo "  $OUT_FASTA"
echo
echo "New BAMs:"
echo "  $OUT_DIR"
echo
echo "Original files were NOT modified."