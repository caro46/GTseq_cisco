# 04_evaluation

Evaluates how well the panel performs after sequencing: read depth, on-target rate, primer dimers, per-locus success, effects of dilution and difference between species. Scripts in `optimization_rounds/` were used to decide which primers to keep or rebalance between rounds. Scripts in `final_panel/` evaluate the finished panel.

## optimization_rounds/

| Script | Round | What it does | Input | Output |
|---|---|---|---|---|
| `gtseq_eval.R` | 1 | Count reads with a single primer or both (on-target) per locus and per sample | FASTQ, primers (primer3 format) | for each sample count table |
| `gtseq_eval_plot.R` | 1 | Plots results of `gtseq_eval.R` to identify bad boys primers, good/bad loci | count tables | figures (.png), list of primers to rebalance or remove |
| `gtseq_eval_sub.sh` | 1 | Submits `gtseq_eval.R` as a cluster job | | |
| `gtseq_eval_optimize.R` | 2-5 | `gtseq_eval.R` taking commandline arguments | FASTQ, primers (IDT format) from previous run (tested and removed) | for each sample count table |
| `gtseq_eval_optimize_plot.R` | 2 | Plots results to identify bad boys, with dilution series | count tables | PDF diagnostic, list of primers to rebalance or remove |
| `gtseq_eval_optimize_plot_no_dil.R` | 3 | Major change = identify multiple class of bad boys (still dilution despite name) | count tables | PDF diagnostic, list of primers to remove in 2 groups for multi testing |
| `gtseq_eval_optimize_sub.sh`, `gtseq_eval_optimize_plot_sub.sh` | 2 | Cluster submission scripts for the two R scripts above | | |
| `gtseq_seq_success_dilution_species.R` | 4 | Sequencing success from multiple primer pools and species | read counts per locus from multiple pool sizes | PDF diagnostic, list of primers to remove |

## final_panel/

| Script | What it does | Input | Output |
|---|---|---|---|
| `gtseq_species_eval_dimer_locus.R` | Per-locus performance and dimer content for each species | read counts per locus | PDF diagnostic |

## Notes

- Round numbers refer to the sequencing rounds of the panel (ex. Round 1 = original pool).
- The three `gtseq_eval*` scripts are variants of the same analysis (depending if multiple pools, dilutions, species were sequenced, which varied between runs).