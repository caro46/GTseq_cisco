# GTseq_cisco

Repo for the GTseq cisco project with Nick Sard. Contains pipelines, scripts, and commands.

Files description:

[To be updated]

## Workflow

| Step | Folder | What it does |
|---|---|---|
| 1 | `scripts/01_discovery/` | Map ascertainment RADseq data, call and filter candidate SNPs |
| 2 | `scripts/02_panel_design/` | Select loci, design primers, prepare orders |
| 3 | `scripts/03_genotyping/` | Extract amplicon reference sequences, map amplicon reads, call genotypes, convert formats |
| 4 | `scripts/04_evaluation/` | Panel performance: depth, missingness, dimers, replicates, polymorphism |

## Software versions

- Mapping and genotyping: bwa-0.7.17, samtools 1.20, gatk4 (v4.2.0.0, HTSJDK Version: 2.24.0, Picard Version: 2.25.0), vcftools 0.1.16, stacks-2.67
- Natural selection: bayescan 2.0
- Primer design: primer3 v. 2.6.1
- R v.4.4.2, Python v.3.9.25 