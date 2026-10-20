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

# ordered primers - to compare names
ordered <- fread("/home2/SardLab/cauret/Ackiss_RADseq_Data/primers_design/ind_filt/relax_criteria/cisco_100wind_relax_primers_ordering_format.tab")

Forward_adapter <- "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
Reverse_adapter <- "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

# for smart_clean_primers.py need:
#- each marker updated F and R primer names
#(if some primers of different markers were the same we only ordered a single one, need to keep good track)


# Make 'Name' the same within each 'Sequence_adapter' group
# by taking the first value within each group
primer_long <- primer %>%
  select(PRIMER_ID,PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE) %>%
  pivot_longer(!PRIMER_ID, names_to = "Primer_orientation", values_to = "Sequence") %>%
    mutate(
    Name = case_when(
      Primer_orientation == "PRIMER_LEFT_SEQUENCE" ~ paste0(PRIMER_ID, "_f1"),
      TRUE ~ paste0(PRIMER_ID, "_r1")
      ),
    Sequence_adapter = case_when(
      Primer_orientation == "PRIMER_LEFT_SEQUENCE" ~ paste0(Forward_adapter,Sequence),
      TRUE ~ paste0(Reverse_adapter,Sequence)
      )
    ) %>%
    group_by(Sequence_adapter) %>%
    mutate(Name = first(Name),
      count=n()) %>% #as a safety check
    ungroup()

head(primer_long)
primer_long %>% filter(!Name %in% ordered$Name) # A tibble: 0 × 6
length(unique(primer_long$Name)) #1266
length(unique(ordered$Name)) #1266
nrow(primer_long)
#[1] 1284
nrow(ordered)
#[1] 1266

primer_wide <- primer_long %>% 
  select(PRIMER_ID, Name, Primer_orientation) %>%
  pivot_wider(
    names_from = Primer_orientation,
    values_from = Name
    )
names(primer_wide) <- gsub("_SEQUENCE", "", names(primer_wide))

# save
write.table(primer_wide, 
  file="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/marker_real_primer_set.txt", 
  row.names=FALSE, col.names=FALSE,
  quote = FALSE, sep="\t") 
