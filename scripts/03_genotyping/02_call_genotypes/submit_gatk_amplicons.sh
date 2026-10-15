#!/bin/bash
# Submit the GATK pipeline to SLURM as a chain of dependent jobs:
#   index -> HaplotypeCaller (array, 1 task per BAM) -> GenomicsDBImport
#         -> GenotypeGVCFs -> hard filtering
# Usage: bash submit_gatk_amplicons.sh
set -euo pipefail

############################ CONFIG (edit me) ############################
export GATK=/programs/gatk4/gatk
export GENOME=/local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons_gatk/cisco_amplicon_gatk.fa
export BAM_DIR=/local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons_gatk            # contains *_sorted.bam
export GVCF_DIR=/local/workdir/cauret/gtseq_cisco/seq_panel_5/GATK_amplicons            # per-sample g.vcf.gz and final VCFs go here
export GDB_PATH=/local/workdir/cauret/gtseq_cisco/seq_panel_5/GATK_amplicons/genomicsdb     # must NOT already exist
export TMP_DIR=/local/workdir/cauret/gtseq_cisco/seq_panel_5/GATK_amplicons/tmp
export INTERVALS=/local/workdir/cauret/gtseq_cisco/seq_panel_5/bam_fq_amplicons_gatk/cisco_amplicon_gatk.GATKintervals.list #Interval list: contig names. [ grep '>' genome.fa | sed 's/>//' | awk '{print $1}' > intervals.list ]
export PREFIX=joint_genotyped           # output name prefix for joint VCF
export NUMB_INT=4                       # intervals imported in parallel (= CPUs for DB step)
LOG_DIR=/local/workdir/cauret/gtseq_cisco/seq_panel_5/GATK_amplicons/slurm_logs
MAX_CONCURRENT=10                       # max HaplotypeCaller tasks running at once

# Optional cluster settings (leave empty if not needed)
PARTITION=""
ACCOUNT=""

# Resources per step (Java heap is set 4-8G below the SLURM memory request)
HC_TIME=04:00:00;     HC_MEM=8G;     HC_CPUS=2;   export HC_HEAP=6G; export HC_CPUS
DB_TIME=04:00:00;     DB_MEM=12G;    DB_CPUS=$NUMB_INT; export DB_HEAP=8G
JOINT_TIME=02:00:00;  JOINT_MEM=12G; JOINT_CPUS=1; export JOINT_HEAP=8G
FILT_TIME=01:00:00;   FILT_MEM=8G;   FILT_CPUS=1; export FILT_HEAP=4G

# Extra HaplotypeCaller arguments. Amplicon data: if your BAMs have duplicates
# MARKED, set this to "--disable-read-filter NotDuplicateReadFilter", otherwise
# HaplotypeCaller discards nearly all reads (amplicon reads look like duplicates).
export HC_EXTRA=""
##########################################################################

mkdir -p "$LOG_DIR" "$GVCF_DIR" "$TMP_DIR"

COMMON=(--parsable --export=ALL --kill-on-invalid-dep=yes --mail-user=caroline.tupper@oswego.edu --mail-type=ALL --nodelist=cbsusard --output="$LOG_DIR/%x_%A_%a.out" --error="$LOG_DIR/%x_%A_%a.err")
[[ -n "$PARTITION" ]] && COMMON+=(--partition="$PARTITION")
[[ -n "$ACCOUNT"   ]] && COMMON+=(--account="$ACCOUNT")

# List of BAMs -> one array task each
export BAM_LIST="$LOG_DIR/bam_list.txt"
ls "$BAM_DIR"/*.bam > "$BAM_LIST"
N_BAM=$(wc -l < "$BAM_LIST")
echo "Found $N_BAM BAM files"
[[ "$N_BAM" -gt 0 ]] || { echo "No *.bam in $BAM_DIR"; exit 1; }

# ---------------- 1. Index reference (only if needed) ----------------
DICT="${GENOME%.*}.dict"
DEP_INDEX=()
if [[ -f "${GENOME}.fai" && -f "$DICT" ]]; then
  echo "Reference index (.fai) and dictionary (.dict) found: skipping indexing"
else
JID_INDEX=$(sbatch "${COMMON[@]}" --job-name=gatk_index --time=01:00:00 --mem=8G --cpus-per-task=1 <<'EOF'
#!/bin/bash
set -euo pipefail
#module load samtools 2>/dev/null || true   # remove/adjust if not needed
[[ -f "${GENOME}.fai" ]] || samtools faidx "$GENOME"
DICT="${GENOME%.*}.dict"
[[ -f "$DICT" ]] || $GATK CreateSequenceDictionary -R "$GENOME"
EOF
)
echo "index:        $JID_INDEX"
DEP_INDEX=(--dependency=afterok:$JID_INDEX)
fi

# ---------------- 2. HaplotypeCaller (array) ----------------
JID_HC=$(sbatch "${COMMON[@]}" ${DEP_INDEX[@]+"${DEP_INDEX[@]}"} --job-name=gatk_hc \
  --array=1-${N_BAM}%${MAX_CONCURRENT} --time=$HC_TIME --mem=$HC_MEM --cpus-per-task=$HC_CPUS <<'EOF'
#!/bin/bash
set -euo pipefail
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
BASE=$(basename "$BAM" .bam)
$GATK --java-options "-Xmx${HC_HEAP}" HaplotypeCaller \
  -R "$GENOME" -I "$BAM" \
  -O "$GVCF_DIR/${BASE}.g.vcf.gz" -ERC GVCF \
  --native-pair-hmm-threads "$HC_CPUS" \
  --max-reads-per-alignment-start 0 \
  $HC_EXTRA \
  --tmp-dir "$TMP_DIR"
EOF
)
echo "HaplotypeCaller: $JID_HC"

# ---------------- 3. GenomicsDBImport ----------------
JID_DB=$(sbatch "${COMMON[@]}" --dependency=afterok:$JID_HC --job-name=gatk_db \
  --time=$DB_TIME --mem=$DB_MEM --cpus-per-task=$DB_CPUS <<'EOF'
#!/bin/bash
set -euo pipefail
V_ARGS=()
for f in "$GVCF_DIR"/*.g.vcf.gz; do V_ARGS+=(-V "$f"); done
$GATK --java-options "-Xmx${DB_HEAP}" GenomicsDBImport \
  --genomicsdb-workspace-path "$GDB_PATH" \
  "${V_ARGS[@]}" \
  --tmp-dir "$TMP_DIR" \
  --max-num-intervals-to-import-in-parallel "$NUMB_INT" \
  --intervals "$INTERVALS"
EOF
)
echo "GenomicsDBImport: $JID_DB"

# ---------------- 4. GenotypeGVCFs ----------------
JID_JOINT=$(sbatch "${COMMON[@]}" --dependency=afterok:$JID_DB --job-name=gatk_joint \
  --time=$JOINT_TIME --mem=$JOINT_MEM --cpus-per-task=$JOINT_CPUS <<'EOF'
#!/bin/bash
set -euo pipefail
$GATK --java-options "-Xmx${JOINT_HEAP}" GenotypeGVCFs \
  -R "$GENOME" -V "gendb://$GDB_PATH" \
  --tmp-dir "$TMP_DIR" \
  -O "$GVCF_DIR/${PREFIX}.vcf.gz"
EOF
)
echo "GenotypeGVCFs:   $JID_JOINT"

# ---------------- 5. Hard filtering ----------------
JID_FILT=$(sbatch "${COMMON[@]}" --dependency=afterok:$JID_JOINT --job-name=gatk_filter \
  --time=$FILT_TIME --mem=$FILT_MEM --cpus-per-task=$FILT_CPUS <<'EOF'
#!/bin/bash
set -euo pipefail
IN="$GVCF_DIR/${PREFIX}.vcf.gz"
OUT="$GVCF_DIR/${PREFIX}"
G="$GATK --java-options -Xmx${FILT_HEAP}"

# ---- SNPs ----
$G SelectVariants -V "$IN" --select-type-to-include SNP \
  -O "${OUT}_snps.vcf.gz"
$G VariantFiltration -V "${OUT}_snps.vcf.gz" \
  -filter "QD < 2.0"  --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "${OUT}_snps_filtered.vcf.gz"
$G SelectVariants --exclude-filtered -V "${OUT}_snps_filtered.vcf.gz" \
  -O "${OUT}_snps_filtered_removed.vcf.gz"
rm -f "${OUT}_snps.vcf.gz"{,.tbi} "${OUT}_snps_filtered.vcf.gz"{,.tbi}

# ---- Indels (MIXED sites are treated with indels) ----
$G SelectVariants -V "$IN" \
  --select-type-to-include INDEL --select-type-to-include MIXED \
  -O "${OUT}_indels.vcf.gz"
$G VariantFiltration -V "${OUT}_indels.vcf.gz" \
  -filter "QD < 2.0"  --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O "${OUT}_indels_filtered.vcf.gz"
$G SelectVariants --exclude-filtered -V "${OUT}_indels_filtered.vcf.gz" \
  -O "${OUT}_indels_filtered_removed.vcf.gz"
rm -f "${OUT}_indels.vcf.gz"{,.tbi} "${OUT}_indels_filtered.vcf.gz"{,.tbi}

# ---- Merge passing SNPs + indels ----

# ---- Sort filtered VCFs before merging ----
$G SortVcf \
  -I "${OUT}_snps_filtered_removed.vcf.gz" \
  -O "${OUT}_snps_filtered_removed.sorted.vcf.gz"

$G SortVcf \
  -I "${OUT}_indels_filtered_removed.vcf.gz" \
  -O "${OUT}_indels_filtered_removed.sorted.vcf.gz"

$G MergeVcfs \
  -I "${OUT}_snps_filtered_removed.sorted.vcf.gz" \
  -I "${OUT}_indels_filtered_removed.sorted.vcf.gz" \
  -O "${OUT}_GATK_hard_filtered_removed.vcf.gz"

rm -f "${OUT}_snps_filtered_removed.sorted.vcf.gz"{,.tbi}
rm -f "${OUT}_indels_filtered_removed.sorted.vcf.gz"{,.tbi}

EOF
)
echo "Hard filtering:  $JID_FILT"
echo "All jobs submitted. Monitor with: squeue -u \$USER"
