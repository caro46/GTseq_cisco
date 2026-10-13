# 02_panel_design

Selects loci from the candidate SNPs, designs and cleans primers, and prepares the order. Includes primer dimer predictions.

## 01_select_loci_make_primers/

| Script | What it does | Input | Output |
|---|---|---|---|
| `vcf_to_primers.R` | Selects SNPs and designs primers around them | candidate-SNP VCF, reference genome | |
| `vcf_to_primers_lib.R`, `vcf_to_primers_fct.R` | Shared functions and libraries sourced by `vcf_to_primers*.R` | | |
| `vcf_to_primers_pop.R` |   | |  |
| `micro_to_primers.R` |  | | primer table |
| `primer3.R` | Fonction to run primer3 in R and parses the results | target sequences | candidate primers |
| `add_single_rep_snp.R` | Adds a single representative SNP for [loci without a multi-SNP window] | | updated locus list |

## 02_primer_cleanup_and_qc/

| Script | What it does | Input | Output |
|---|---|---|---|
| `clean_primers.py`, `smart_clean_primers_v1.py` | Removes primers that fail delta G filter (proxy for dimers) | mfeprimer or mfeprimer/primerpooler outputs | cleaned primer table |
| `primers_tab_to_fasta.R` | Converts the primer table to FASTA | primer table | FASTA |
| `primer_tab_to_real_list.R` | Builds the final primer list for ordering | cleaned primer table | primer list |
| `find_share_primers.R` | Finds primers shared between studies | primer table | table of shared primers |
| `compare_mfe_pooler.R` | Compares minimum free energy / pooling predictions | mfeprimer, primerpooler outputs | comparison figs |
| `parse_dimer_report.py` | Parses the primer dimer report (mfeprimer, primerpooler)| dimer report | CSV |
| `raw_dimer_report_compare.R` | Compares dimer reports between software | parsed dimer tables | comparison, plots |


## Software

R v.4.4.2, Python v.3.9.25; primer3 v. 2.6.1; MFEprimer v.3.2.6, PrimerPooler v.1.84

## Notes

### Commands for mfeprimer and primerpooler:

`mfeprimer dimer -i ../../cisco_100wind_relax_primers.fa -o cisco_100wind_relax_primers_dimer_GTseq_param_mfe_tm57.txt --mono 50 --diva 3 --dntp 0.2 --oligo 54 --tm 57`

`pooler --pools=2,,poolerfile /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_primerpooler.fa --dg=330.15,3,50,0.2 --print-bonds=-2 >pooler_bonds.txt`
