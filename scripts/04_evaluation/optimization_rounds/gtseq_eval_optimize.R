#module load R/4.4.2
packages <- c("Biostrings","dplyr", "data.table", "ggplot2", "tidyr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")
library(ShortRead)
# Packages loading
lapply(packages, library, character.only = TRUE)

# Access the Task ID 
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

# Other command line arguments
## Primer in IDT ordering format
primer_table_file <- args[2]
## Path to fastq
fastq_dir <- args[3]
## Output directory
out_eval <- args[4]
## Primers to remove
to_remove_file <- args[5]
## Max number of mismatched allowed in primer regions
p_mismatch = as.numeric(args[6])
##Suffix to define forward and reverse primers
fwd_suffix <- "_[fF]"
rev_suffix <- "_[rR]"

# List all fastq.gz files
#fastq_files <- list.files(path = fastq_dir, pattern = "\\R1.fastq\\.gz$", full.names = TRUE)
fastq_files <- list.files(path = fastq_dir, pattern = "\\R1.fastq(\\.gz)?$", full.names = TRUE)

## Loading inputs
primer_table <- fread(primer_table_file)
to_remove <- fread(to_remove_file)

# Prepare primer table in good format
primers_tested <- primer_table %>% filter(! PRIMER_ID %in% to_remove$primer_set)

## Test
# primer_dir <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/"
# output_filename_keep <- paste0(primer_dir,"primers_to_keep_for_run2.csv")
# #primer_table <- fread(output_filename_keep)
# primers_run1 <- fread("/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/tm57/cisco_100wind_relax_primers_ordering_format_all_dimers.csv", header = TRUE, sep = ",")
# output_filename_remove <- paste0(primer_dir,"primers_to_remove_for_run2.csv")
# to_remove <- fread(output_filename_remove)
# primers_run1_tested <- primers_run1 %>% filter(! PRIMER_ID %in% to_remove$primer_set)
# out_eval <- "/local/workdir/cauret/gtseq_cisco/seq_panel_1/eval/"
# p_mismatch = 3
# primers_tested <- primers_run1_tested

#col to have: PRIMER_ID, PRIMER_LEFT_SEQUENCE, PRIMER_RIGHT_SEQUENCE

#case_when to add column, _f = PRIMER_LEFT_SEQUENCE; _r = PRIMER_RIGHT_SEQUENCE
primers_tested <- primers_tested %>%
  mutate(
    primer_direction = case_when(
      grepl(fwd_suffix, Name) ~ "PRIMER_LEFT_SEQUENCE", 
      grepl(rev_suffix, Name) ~ "PRIMER_RIGHT_SEQUENCE",
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
  )

# Safety check: detect incomplete pairs
primers_tested_wide <- primers_tested_wide %>%
filter(
  ! is.na(PRIMER_LEFT_SEQUENCE) |
  ! is.na(PRIMER_RIGHT_SEQUENCE)
  )


# Function to check primers in one file for one primer set
check_primers <- function(file, set_name, fwd_seq, rev_seq) {
  fq <- readFastq(file)
  seqs <- sread(fq)
  
  fwd <- DNAString(fwd_seq)
  rev <- DNAString(rev_seq)
  fwd_rc <- reverseComplement(fwd)
  rev_rc <- reverseComplement(rev)

  has_forward <- vcountPattern(fwd, seqs, max.mismatch=p_mismatch) > 0
  has_reverse <- vcountPattern(rev, seqs, max.mismatch=p_mismatch) > 0
  has_forward_rc <- vcountPattern(fwd_rc, seqs, max.mismatch=p_mismatch) > 0
  has_reverse_rc <- vcountPattern(rev_rc, seqs, max.mismatch=p_mismatch) > 0

  
  # to add: fwd_rc and rev; fwd and rev_rc; fwd_rc and rev_rc; remove neither
  #only_forward <- has_forward & !has_reverse
  #only_reverse <- has_reverse & !has_forward
  #both_primers <- has_forward & has_reverse
  #neither <- !has_forward & !has_reverse
  only_forward <- has_forward & !has_reverse & !has_reverse_rc & !has_forward_rc
  only_reverse <- has_reverse & !has_forward & !has_reverse_rc & !has_forward_rc
  only_forward_rc <- has_forward_rc & !has_reverse & !has_reverse_rc & !has_forward
  only_reverse_rc <- has_reverse_rc & !has_forward & !has_reverse & !has_forward_rc
  #both_primers <- has_forward & has_reverse
  forward_reverse <- has_forward & has_reverse
  forward_rc_reverse <- has_forward_rc & has_reverse
  forward_reverse_rc <- has_forward & has_reverse_rc
  forward_rc_reverse_rc <- has_forward_rc & has_reverse_rc
  #neither <- !has_forward & !has_reverse
  
  
  # data.frame(
  #   file = basename(file),
  #   primer_set = set_name,
  #   total_reads = length(seqs),
  #   only_forward = sum(only_forward),
  #   only_reverse = sum(only_reverse),
  #   both_primers = sum(both_primers),
  #   neither = sum(neither)
  # )
    data.frame(
    file = basename(file),
    primer_set = set_name,
    total_reads = length(seqs),
    only_forward = sum(only_forward),
    only_reverse = sum(only_reverse),
    only_forward_rc = sum(only_forward_rc),
    only_reverse_rc = sum(only_reverse_rc),
    forward_reverse = sum(forward_reverse),
    forward_rc_reverse = sum(forward_rc_reverse),
    forward_reverse_rc = sum(forward_reverse_rc),
    forward_rc_reverse_rc = sum(forward_rc_reverse_rc)
    #neither = sum(neither)
  )
}

# Loop through all combinations of files and primer sets
results <- list()
#for (file in fastq_files) {
#for (file in fastq_files[1:5]){ #for now subset otherwise taking forever and not necessary to see if primers ok
for (file in fastq_files[task_id]){
  #for (i in 1:2) {
  #for (i in seq_len(nrow(primer_table))) {
    #row <- primer_table[i, ]
  for (i in seq_len(nrow(primers_tested_wide))) {
    row <- primers_tested_wide[i, ]
    #res <- check_primers(file, row$set_name, row$forward, row$reverse)
    res <- check_primers(file, row$PRIMER_ID, row$PRIMER_LEFT_SEQUENCE, row$PRIMER_RIGHT_SEQUENCE)
    results <- append(results, list(res))
  }
}


# Combine all results into dataframe
results_df <- bind_rows(results) %>%
  mutate(total_primers = rowSums(across(only_forward:forward_rc_reverse_rc)),
        individuals = gsub(".fastq.gz", "", file))

# Create a filename using the task ID
output_filename <- paste0(out_eval,"primers_reads_count_", task_id, ".txt")

# Write data from one individual to a file
write.csv(results_df, file = output_filename, row.names = FALSE, quote = FALSE)

