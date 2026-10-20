#module load R/4.4.2
packages <- c("Biostrings","dplyr", "data.table", "ggplot2")

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


# Prepare primer table in good format
out_wind <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/primers_design/ind_filt/relax_criteria/"
primer_temp <- fread(paste(out_wind,"cisco_100wind_relax_master_filtered_final.tab",sep=""))
to_remove <- fread(paste(out_wind,"remove.these.uHapV1.markers.txt",sep="")) %>%
  dplyr::rename("PRIMER_ID" = "V1", "Ampl_numb" = "V2")
primer_table <- primer_temp %>% dplyr::filter(!PRIMER_ID %in% to_remove$PRIMER_ID)
#cisco_100wind_relax_master_filtered_final_no_offtarget.tab

# Load primer table
#primer_table <- read.csv("path/to/primer_table.csv", stringsAsFactors = FALSE)

# List all fastq.gz files
#fastq_files <- list.files(path = "/local/workdir/cauret/gtseq_cisco/seq_panel_1/temp", pattern = "\\.fastq\\.gz$", full.names = TRUE)
fastq_files <- list.files(path = "/local/workdir/cauret/gtseq_cisco/seq_panel_1", pattern = "\\.fastq\\.gz$", full.names = TRUE)

#fastq_files <- list.files(path = "path/to/fastq/files", pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Output directory
out_eval <- "/local/workdir/cauret/gtseq_cisco/seq_panel_1/eval/"

# variables to specify
## Max number of mismatched allowed in primer regions
p_mismatch = 3

# need to add $SLURM_ARRAY_TASK_ID number as a numeric for chosing file

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
  for (i in seq_len(nrow(primer_table))) {
    row <- primer_table[i, ]
    #res <- check_primers(file, row$set_name, row$forward, row$reverse)
    res <- check_primers(file, row$PRIMER_ID, row$PRIMER_LEFT_SEQUENCE, row$PRIMER_RIGHT_SEQUENCE)
    results <- append(results, list(res))
  }
}


# Combine all results into dataframe
results_df <- bind_rows(results) %>%
  mutate(total_primers = rowSums(across(only_forward:forward_rc_reverse_rc)),
        individuals = gsub(".fastq.gz", "", file))
#print(results_df)
# write.table(results_df, 
#   #file="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_reads_count_5ind.txt", row.names=FALSE,
#   file="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_reads_count_.txt", row.names=FALSE,
#   quote = FALSE) 

# Create a filename using the task ID
output_filename <- paste0(out_eval,"primers_reads_count_", task_id, ".txt")

# Write data from one individual to a file
write.csv(results_df, file = output_filename, row.names = FALSE, quote = FALSE)



# heatmap with markers at bottom and ind on side, color = numb reads
#maybe add colm that has any combination of forw and rev?
# need to remove file ext from ind, change orientation of name on x axis 
# hp_read_count <- ggplot(results_df, 
#   aes(x = primer_set, y = individuals, fill = total_primers,color = "black")) +
#   geom_tile(color = "black",
#             #color = "white",
#             #lwd = 1.5,
#             linetype = 1
#             ) +
#   scale_fill_gradient(low = "white", high = "red") +
#   #coord_fixed() +
#   theme(axis.text.x = element_blank()#, # Remove x-axis tick labels
#         #axis.text.y = element_blank()
#         )#,
#         #axis.title.x = element_blank(), # Remove x-axis title
#         #axis.title.y = element_blank()) # Remove y-axis title # Remove y-axis tick labels
# #ggsave(hp_read_count, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/hp_read_count.png")
# ggsave(hp_read_count,
#   width = 20,
#   height = 6,
#   units = "in", 
#   filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/hp_read_count.png")

#checking issue (most only f or rev)

# primer_table <- primer_temp %>% dplyr::filter(PRIMER_ID == "60028_63")
# # previously only_forward 25 reads, will check on those reads
# fq <- readFastq("/local/workdir/cauret/gtseq_cisco/seq_panel_1/temp/CT_cART_M2Q_10483695_M53J6_R1.fastq.gz")
# seqs <- sread(fq)
# fwd <- DNAString(primer_table$PRIMER_LEFT_SEQUENCE)
# rev_rc <- reverseComplement(DNAString(primer_table$PRIMER_RIGHT_SEQUENCE))
# has_forward <- vcountPattern(fwd, seqs, max.mismatch=p_mismatch) > 0
# # Identify sequences where the pattern occurs at least once (count > 0)
# sequences_with_pattern <- seqs[has_forward]
# as.character(sequences_with_pattern[1])
# #[1] "CCAAGATTCTTTGCACTCTG GGAGGAGGACACAATGGAGTTGTCAACCGAGATGG tail> AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAATAAATAAAATAAAAGTAGTGAGTTTGAAGAATGGTAGATAAATTATAGGGAAGAAAATGTGAAATAGTTAATAGACGGAGAGGATTGTAGACTATTGTTTGGTGAAGAGTGTGTATTGGTGTAGTTGAATGTTGAGGTGTAGTATTGTAGGTCAGTATGATGTATGTGTCGATTGT"
# #> primer_table
# #   PRIMER_ID PRIMER_NUM   PRIMER_LEFT_SEQUENCE PRIMER_RIGHT_SEQUENCE
# #      <char>      <int>                 <char>                <char>
# #1:  60028_63          1 CCAAGGTTCTTTGTACTCTGGG    AGCTCCTCTTCCTCCCGG

# #From ordering table
# #2114753_63  PRIMER_RIGHT_SEQUENCE CGATCTCTAGCTGCTGATGGTAC 2114753_63_r1 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGATCTCTAGCTGCTGATGGTAC 1

# #as.character(sequences_with_pattern[2])
# #[1] "CCAAGGTTCTTTGCACTCTG GGAGGAGGACACAATGGAGTTGTCAACCGAGATGG tail> AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTCGATCTCGTATGCCGTCTTCTGCTTTAAAAAAAAAAAAAAAAAAAAAAATATAAATTAGGTATGTAGATGAGGAGATAAAGGAGAAAAAATAGTGAGTTAGGAAGATGTAGTGATTAGTTTGGAGTGGGAGTGTTAGGGATGTTAGGAAGAGGTAGGTGTAGATTGAAAGAGTAGGAGAGTATATTCCATGTATAGTGTGGAGAGATATAGGTTGTT"

# # Get the length of each read
# read_lengths <- width(seqs)
# length_counts <- table(read_lengths)
# print(length_counts)

# #same with 35507_18
# primer_table <- primer_temp %>% dplyr::filter(PRIMER_ID == "35507_18")
# seqs <- sread(fq)
# fwd <- DNAString(primer_table$PRIMER_LEFT_SEQUENCE)
# rev_rc <- reverseComplement(DNAString(primer_table$PRIMER_RIGHT_SEQUENCE)) #GGCGGGTCATGTTTTGGT
# has_reverce_rc <- vcountPattern(rev_rc, seqs, max.mismatch=p_mismatch) > 0
# # Identify sequences where the pattern occurs at least once (count > 0)
# sequences_with_pattern <- seqs[has_reverce_rc]
# as.character(sequences_with_pattern[1])
# #[1] "TT rev_rc GGTGGGTCATGTTTCGGA GG tail>AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC GGCTCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAATAAACAAATATATGAAAATAAGAGTAAAGAAAAAAATAAGAGGAGAAGAAGAGATGAAGAGAGATGAAGTAGAGTTATAGGTAAGGATTTATTAAAATAGGTATAGTCTAATATTATAATTGGGGTGGGTTTAGTGTATGATGAGTTTTTTTTGTATATGAATGAGATGATATGGGAGATAGTATGAGTATAGTAGTAGATGATATAGTT"

#  as.character(sequences_with_pattern[3])
# #[1] "TTGGTGGGTCATGTTTCGGG AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC GGCTCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATAAAAAAAAGAAAAAGGAAGGAAAAAAAAAAAAGAGAAAATAAAAAGAAAGAGAAGATAAGAGAAAGTGGGGAGGAGAGAAGAAGGAAGAAGGGGAAGAATAGGTATAAGGAGGTAAAGAGAAAGAAAATGTGAGAGATATAATGAGAGTGAGTAAAGTTAAGAAATAAGAAGTAGGTAAGAGAGTATATAGATAGTGAAAGTAAGGA"

# #cRead1: 5'-CGACAGGTTCAGAGTTCTACAGTCCGACGATC-3'
# #cIndex2: 5’-GATCGTCGGACTGTAGAACTCTGAACCTGTCG-3’
# #9:38
# #cIndex1: 5'-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3’
# #cRead2: 5'-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3'

# matches <- matchPattern(pattern = rev_rc, subject = as.character(sequences_with_pattern[1]),max.mismatch=3)

