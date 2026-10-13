# 01_discovery

Uses the RAD-seq data to find candidate SNPs for the panel. The main output is a filtered VCF of candidate SNPs, which is the input for `02_panel_design/`.


| Script | What it does | Input | Output |
|---|---|---|---|
| `fasta_to_vcf.py` | Maps paired-end reads with bwa, sorts and indexes with samtools, calls genotypes with GATK or stacks | FASTQ, reference genome | sorted BAM, VCF |
| `extract_bam.py` | Extracts reads or regions from BAM files | BAM, [region list] | subset BAM |
| `extract_bam_sub.sh` | Cluster submission for `extract_bam.py` | | |
| `vcf_explo.R` | Summaries of the VCF: depth, missingness, ...| VCF | tables, plots |


## Software

- bwa-0.7.17, samtools 1.20, gatk4 (v4.2.0.0, HTSJDK Version: 2.24.0, Picard Version: 2.25.0), vcftools 0.1.16, stacks-2.67
- R v.4.4.2, Python v.3.9.25

## Notes
