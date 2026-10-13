# 03_genotyping

Turns amplicon sequencing reads into genotype calls. Converts into SNP matrix, compatible with Cervus and Sequoia to run parentage analysis.

## Scripts

| Script | What it does | Input | Output |
|---|---|---|---|
| `01_trim_and_map/trim_by_primers.py` | Trim FASTA sequences using primer pairs and generate genomic coordinate annotations | FASTA catalog, Primers | FASTA amplicon seq, genomic coordinate table |
| `01_trim_and_map/bwa_amplicon.sh` | Maps reads to the amplicon reference and extract mapped reads (targets) | amplicon FASTQ, amplicon reference sequences | target FASTQ, BAM |
| `02_call_genotypes/rename_files_issues_GATK.sh` | Rename the amplicon sequence compatible with GATK | FASTA, BAM | FASTA, BAM |
| `02_call_genotypes/submit_gatk_amplicons.sh` | Calls variants from amplicon data with GATK | BAM | VCF |
| `02_call_genotypes/gtseq_microhap*_sub.sh.sh` | Calls variants from amplicon data with gtseq_microhap | FASTQ, primers | VCF |
| `03_convert_formats/microhap_vcf_to_database.R` | Converts VCF to genotype matrix | VCF | SNP matrix (.csv) |

## Software and settings

bwa-0.7.17, samtools 1.20, gatk4 (v4.2.0.0, HTSJDK Version: 2.24.0, Picard Version: 2.25.0). Key parameters: GATK `--max-reads-per-alignment-start 0` (disables downsampling for reads that share the same alignment start position).

## Notes

Scripts used on last panel round (optimized panel). The output from `trim_by_primers.py` is not directly compatible for `GATK`, use `rename_files_issues_GATK.sh` to make it compatible. However because of that, note that the locus naming between `GATK` and `gtseq_microhap` might not be exactly the same (ex. 1004433_1 = original catalog seq name `_` SNP number; 1004433_23_1 = locus name in primer files or original catalog seq name `_` window start `_` SNP number). 