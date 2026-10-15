#module load R/4.4.2
packages <- c("dplyr", "data.table", "tidyr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)

# Access the Task ID 
args <- commandArgs(trailingOnly = TRUE)

# Other command line arguments
## Primer in IDT ordering format
primer_table_file <- args[1]
## Primers to remove
to_remove_file <- args[2]
## Primer in locus format
locus_primer <- args[3]
## Output directory
out_dir <- args[4]
## Outfile name
out_file <- args[5]
# Create a filename
output_filename <- paste0(out_dir,out_file)


## Test
primer_table <- fread("/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/tm57/cisco_100wind_relax_primers_ordering_format_all_dimers.csv", header = TRUE, sep = ",")
to_remove_file ="/local/workdir/cauret/gtseq_cisco/seq_panel_4/primers_to_remove_combined_loci_367.csv"
to_remove <- fread(to_remove_file)
output_filename <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_5_files/locus_primers_367.csv"

##Suffix to define forward and reverse primers
fwd_suffix <- "_[fF]"
rev_suffix <- "_[rR]"

## Loading inputs
primer_table <- fread(primer_table_file)
to_remove <- fread(to_remove_file)

# Prepare primer table in good format
primers_tested <- primer_table %>% filter(! PRIMER_ID %in% to_remove$primer_set)

primers_tested <- primers_tested %>%
  mutate(
    primer_direction = case_when(
      grepl(fwd_suffix, Name) ~ "fwd_primer", 
      grepl(rev_suffix, Name) ~ "rev_primer",
      #grepl("_f", Name) ~ "PRIMER_LEFT_SEQUENCE", 
      #grepl("_r", Name) ~ "PRIMER_RIGHT_SEQUENCE", 
      TRUE ~ "NA"
      )
    )

# Convert to wide format
primers_tested_wide <- primers_tested %>%
  select(PRIMER_ID, primer_direction, Sequence) %>%
  pivot_wider(
    id_cols   = PRIMER_ID,
    names_from = primer_direction, # Column whose values will become new column names
    values_from = Sequence # Column containing the values for the new columns
  ) %>%
  dplyr::rename(locus = PRIMER_ID)


# Write primers to a file
write.csv(primers_tested_wide, file = output_filename, row.names = FALSE, quote = FALSE)


