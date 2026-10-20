#module load R/4.4.2
packages <- c("dplyr", "data.table","tidyr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)

# Prepare primer table in good format
out_wind <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/primers_design/ind_filt/relax_criteria/"
primer <- fread(paste(out_wind,"cisco_100wind_relax_master_filtered_final.tab",sep=""))
to_remove <- fread(paste(out_wind,"remove.these.uHapV1.markers.txt",sep="")) %>%
  dplyr::rename("PRIMER_ID" = "V1", "Ampl_numb" = "V2")
primer <- primer %>% dplyr::filter(!PRIMER_ID %in% to_remove$PRIMER_ID) 

# for primerpooler need:
#- primers seq in fasta format
#- each couple F/R needs to be together
#- Names of amplicons’ primers should end with F or R (for Forward and Reverse), and otherwise match

#primer <- primer %>%
primer_fa <- primer %>%
  select(PRIMER_ID,PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE) %>%
  pivot_longer(!PRIMER_ID, names_to = "Primer_orientation", values_to = "Sequence") %>%
    mutate(
    Name = case_when(
      Primer_orientation == "PRIMER_LEFT_SEQUENCE" ~ paste0(PRIMER_ID, "_F"),
      TRUE ~ paste0(PRIMER_ID, "_R")
      )
    ) %>%
    mutate(Name = paste0(">", Name)) %>%
    select(PRIMER_ID,Name,Sequence) %>%
    pivot_longer(!PRIMER_ID, names_to = "identifier", values_to = "Header_Sequence") %>%
    select(Header_Sequence)

# save
#write.table(primer,
write.table(primer_fa, 
  file="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_primerpooler.fa", 
  row.names=FALSE, col.names=FALSE,
  quote = FALSE) 

# to process output from primerpooler with smart_clean_primers.py needs 
# a tab file with locus_name \t primer1 \t primer2 (=PRIMER_ID with _F or _R added)
primer_tab <- primer %>%
  select(PRIMER_ID) %>%
  mutate(
    Primer_F = paste0(PRIMER_ID, "_F"),
    Primer_R = paste0(PRIMER_ID, "_R"))

write.table(primer_tab, 
  file="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_primerpooler.tab", 
  row.names=FALSE, col.names=FALSE,
  sep = "\t", 
  quote = FALSE) 

