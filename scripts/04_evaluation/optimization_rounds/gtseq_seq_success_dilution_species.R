
### To do:
# add colors on single primer plot for the over-represented 

#module load R/4.4.2
packages <- c("dplyr", "data.table", "ggplot2","stringr","VennDiagram")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
	install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)

# Access Command Arguments 
args <- commandArgs(trailingOnly = TRUE)
out_eval <- args[1]
out_dir <- args[2] #output dir for diag plots and primer list to remove and to keep
cum_limit <- as.numeric(args[3])
run_number <- as.numeric(args[4])

# Output file names
out_dir <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_4_files/"

output_filename_diag_plot <- paste0(out_dir,"primers_panels_seq_diagnostic.pdf")
output_filename_keep <- paste0(out_dir,"primers_to_keep_for_run",run_number,".csv")
#output_filename_remove <- paste0(out_dir,"primers_to_remove_for_run",run_number,".csv")
output_filename_remove <- paste0(out_dir,"primers_to_remove_from_pool399",".csv")
#output_filename_worst <- paste0(out_dir,"primers_worst_for_run",run_number,".csv")

#paste0(out_dir,"primers_to_keep_for_run2.csv")

# Output directorys
out_eval <- c(
	"/local/workdir/cauret/gtseq_cisco/seq_panel_4/eval",
	"/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/299loci/eval",
	"/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/349loci/eval",
	"/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/299loci/eval",
	"/local/workdir/cauret/gtseq_cisco/seq_panel_4/link2/399loci/eval"
	)

# List all read count files
#count_files <- list.files(path = out_eval, pattern = "\\.txt$", full.names = TRUE)

count_files <- unlist(lapply(out_eval, function(dir) {
	list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)
}))

# Read each CSV file into a list of data frames
list_of_dfs <- lapply(count_files, read.csv)

# Combine all data frames into a single data frame
count_table <- do.call(rbind, list_of_dfs) %>%
	select(-file) %>%
	#add primer_set column
	#mutate(loci_number_pool = str_split(individuals, "_", simplify = TRUE)[,2])
	mutate(
		loci_number_pool = str_extract(individuals, "CT\\d+"),
		loci_number_pool = str_replace(loci_number_pool, "CT499", "CT399"),
		loci_number_pool = str_remove(loci_number_pool, "CT"),
		loci_number_pool = as.numeric(loci_number_pool)
) %>%
  mutate(
	rc_bin = case_when(
	  forward_reverse_rc == 0 ~ "0",
	  forward_reverse_rc >= 1 & forward_reverse_rc <= 9 ~ "1-9",
	  forward_reverse_rc >= 10 & forward_reverse_rc <= 99 ~ "10-99",
	  forward_reverse_rc >= 100 & forward_reverse_rc <= 999 ~ "100-999",
	  forward_reverse_rc >= 1000 & forward_reverse_rc <= 9999 ~ "1000-9999",
	  forward_reverse_rc >= 10000 ~ "10000+",
	  TRUE ~ NA_character_
	),
	rc_bin = factor(
	  rc_bin,
	  levels = c("0", "1-9", "10-99", "100-999", "1000-9999", "10000+")
	)
  )

# View the first few rows of the combined table
head(count_table)

## Looking at dimers across pools

#limit_read = 1000000
single_count_plot_fct <- function(data_frame_single_primer) {
  ggplot(data_frame_single_primer, aes(x = primer_set, y = total_primers_single)) +
  geom_point() +
  labs(title = "Primer Set - Single Primer Reads Count (Individuals Sum)",
	   x = "Primer Sets",
	   y = "Total Single Primers") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
		axis.ticks.x = element_blank()) 
}

sum_single_pool <- count_table %>%
	filter(! str_detect(individuals, "_(1x|2x|4x)_")) %>%
	group_by(primer_set, loci_number_pool) %>%
	summarise(
		total_primers_single = sum(
			c_across(where(is.numeric) & contains("only_")),
			na.rm = TRUE
			),
			.groups = "drop"   # 👈 IMPORTANT
			) %>%
	arrange(desc(total_primers_single)) %>%
	mutate(
		primer_set = factor(primer_set, levels = unique(primer_set))
		)

single_count_plot_pool <- single_count_plot_fct(sum_single_pool) + #facet_wrap(~dilution_category)
	facet_grid(
	loci_number_pool ~ .,   # <-- stacked on top of each other
	scales = "free_y",
	space  = "free_y"
	) 

# Subset to focus on dilution effect

dilution_subset <- count_table %>%
	filter(str_detect(individuals, "_(1x|2x|4x)_")) %>%
	rename(files_no_ext = individuals) %>% 
	mutate(
		dilution_factor = str_extract(files_no_ext, "(?<=_)\\d+x(?=_)"),
		# Pattern pieces
		# \\d+ → one or more digits (1, 2, 10, etc.)
		# x → literal “x”
		# (?<=_) → must be preceded by _
		# (?=_) → must be followed by _
		dilution_factor = factor(dilution_factor, levels = c("1x", "2x", "4x")),
		species = case_when(
			str_detect(files_no_ext, "cHOY") ~ "cHOY",
			str_detect(files_no_ext, "cART") ~ "cART",
			TRUE ~ NA_character_
		),
		# NEW individuals column (name before dilution factor)
		individuals = str_match(files_no_ext, "CT\\d+_(.*?)_(1x|2x|4x)")[,2],    
		# Extract cell number
		cell_number = str_extract(files_no_ext, "^\\d+"),
		cell_number = as.numeric(cell_number)
	)

head(dilution_subset)

# Success within a sample
success_by_sample <- dilution_subset %>%
	mutate(success = forward_reverse_rc >= 1) %>%
	group_by(
		individuals,
		cell_number,
		dilution_factor,
		loci_number_pool,
		species
	) %>%
	summarise(
		successful_loci_n = sum(success, na.rm = TRUE),
		total_loci = n(),
		prop_successful = successful_loci_n / total_loci,
		.groups = "drop"
	)

# Overall effect
boxplot_overall <- ggplot(success_by_sample,
			 aes(x = dilution_factor, y = successful_loci_n)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(aes(color = species),
							width = 0.15,
							alpha = 0.7,
							size = 2) +
	facet_grid(species ~ loci_number_pool) +
	labs(
		x = "Dilution factor",
		y = "Number of successful loci (>=1 read)",
		title = "Effect of dilution on number of successfully sequenced loci (Pool 349)"
	) +
	theme_bw()

# Within-individual effect
paired_line_plot <- ggplot(success_by_sample,
			 aes(x = dilution_factor,
					 y = successful_loci_n,
					 group = individuals)) +
	geom_line(alpha = 0.4) +
	geom_point(aes(color = species), size = 2) +
	facet_grid(species ~ loci_number_pool) +
	labs(
		x = "Dilution factor",
		y = "Number of successful loci (>=1 read)",
		title = "Within-individual change in successful loci across dilutions (Pool 349)"
	) +
	theme_bw()

# Focusing on 1x

heatmap_df <- dilution_subset %>%
	filter(dilution_factor == "1x") %>%
	mutate(success = forward_reverse_rc >= 1)

## Order loci by overall success

heatmap_df <- heatmap_df %>%
	group_by(primer_set) %>%
	mutate(locus_order = mean(success)) %>%
	ungroup()

heatmap_df$primer_set <- reorder(heatmap_df$primer_set, heatmap_df$locus_order)

## Order individuals by success rate

heatmap_df <- heatmap_df %>%
	group_by(individuals) %>%
	mutate(ind_order = mean(success)) %>%
	ungroup()

heatmap_df$individuals <- reorder(heatmap_df$individuals, heatmap_df$ind_order)

heatmap_plot <- ggplot(heatmap_df,
			 aes(x = individuals, y = primer_set, fill = forward_reverse_rc)) +
	geom_tile() +
	#scale_fill_viridis_c(trans = "log1p") +
	scale_fill_gradient(
		low = "white",
		high = "red",
		na.value = "grey90"
	) +
	facet_wrap(~ species, scales = "free_x") +
	labs(
		x = "Individual",
		y = "Locus",
		fill = "Read count",
		title = "Locus-level read counts (1x dilution, pool 349)"
	) +
	theme_bw() +
	theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
		axis.text.y = element_blank()#,
		#axis.text.y = element_text(size = 6)
	)


#Code: create range bins
heatmap_df <- dilution_subset %>%
  filter(dilution_factor == "1x") 

#Order loci by mean read count
locus_levels <- heatmap_df %>%
  group_by(primer_set) %>%
  summarise(mean_rc = mean(forward_reverse_rc, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rc) %>%
  pull(primer_set)

#Order individuals by mean read count
individual_levels <- heatmap_df %>%
  group_by(individuals) %>%
  summarise(mean_rc = mean(forward_reverse_rc, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rc) %>%
  pull(individuals)

#Apply ordering
heatmap_df <- heatmap_df %>%
  mutate(
	primer_set = factor(primer_set, levels = locus_levels),
	individuals = factor(individuals, levels = individual_levels)
  )

#Heatmap code
heatmap_plot <- ggplot(heatmap_df,
	   aes(x = individuals, y = primer_set, fill = rc_bin)) +
  geom_tile() +
  scale_fill_manual(
	values = c(
	  "0" = "white",
	  "1-9" = "#FEE5D9",
	  "10-99" = "#FCBBA1",
	  "100-999" = "#FC9272",
	  "1000-9999" = "#DE2D26",
	  "10000+" = "#67000D"
	),
	na.value = "grey90"
  ) +
  facet_wrap(~ species, scales = "free_x") +
  labs(
	x = "Individual",
	y = NULL,
	fill = "Read count",
	title = "Locus-level read count categories at 1x dilution (Pool 349)"
  ) +
  theme_bw() +
  theme(
	axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
	axis.text.y = element_blank(),
	axis.ticks.y = element_blank()
  )

# Now all other individuals
dilution_excluded <- count_table %>%
	filter(! str_detect(individuals, "_(1x|2x|4x)_")) 

#Order loci by mean read count
locus_levels <- dilution_excluded %>%
  group_by(primer_set) %>%
  summarise(mean_rc = mean(forward_reverse_rc, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rc) %>%
  pull(primer_set)

#Order individuals by mean read count
individual_levels <- dilution_excluded %>%
  group_by(individuals) %>%
  summarise(mean_rc = mean(forward_reverse_rc, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rc) %>%
  pull(individuals)

#Apply ordering
dilution_excluded <- dilution_excluded %>%
  mutate(
	primer_set = factor(primer_set, levels = locus_levels),
	individuals = factor(individuals, levels = individual_levels)
  )

heatmap_plot_no_dil <- ggplot(dilution_excluded,
			 #aes(x = individuals, y = primer_set, fill = forward_reverse_rc)) +
			 aes(x = individuals, y = primer_set, fill = rc_bin)) +
	geom_tile() +
	#scale_fill_viridis_c(trans = "log1p") +
	# scale_fill_gradient(
	# 	low = "white",
	# 	high = "red",
	# 	na.value = "grey90"
	# ) +
	scale_fill_manual(
	values = c(
	  "0" = "white",
	  "1-9" = "#FEE5D9",
	  "10-99" = "#FCBBA1",
	  "100-999" = "#FC9272",
	  "1000-9999" = "#DE2D26",
	  "10000+" = "#67000D"
	),
	na.value = "grey90"
  ) +
	facet_wrap(~ loci_number_pool, scales = "free_x") +
	labs(
		x = "Individual",
		y = "Locus",
		fill = "Read count",
		title = "Locus-level read counts"
	) +
	theme_bw() +
	theme(
		#axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
		axis.text.y = element_blank(),
		axis.text.x = element_blank()
		#axis.text.y = element_text(size = 6)
	)

### Add color on Single primer plot based on over-represented both primers
#list of over-represented FR primers
FR_399_high_list <-  unique(as.character( (dilution_excluded %>% 
	filter(forward_reverse_rc >= 1000 & loci_number_pool == 399) )$primer_set))

sum_single_pool <- sum_single_pool %>% mutate(FR_399_high =  case_when(
		primer_set %in% FR_399_high_list ~ "yes",
		TRUE ~ "no"))

single_count_plot_pool_2 <- ggplot(sum_single_pool, 
	aes(x = primer_set, y = total_primers_single, color = FR_399_high)) +
  geom_point() +
  scale_color_manual(
  	values = c(
  		"no"  = "black",
  		"yes" = "purple"
  		)
  	) +
  labs(title = "Primer Set - Single Primer Reads Count (Individuals Sum)",
	   x = "Primer Sets",
	   y = "Total Single Primers") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
		axis.ticks.x = element_blank()) + 
	facet_grid(
	loci_number_pool ~ .,   # <-- stacked on top of each other
	scales = "free_y",
	space  = "free_y"
	) +   
	geom_hline(yintercept = 50000, col = "red", linetype = "dashed", size = 1)


single_399_high_list <-  unique(as.character( (sum_single_pool %>% 
	filter(total_primers_single > 50000 & loci_number_pool == 399) )$primer_set))

merged_high_list <- unique( c(FR_399_high_list, single_399_high_list) )

# Export file with explicit primer names
merged_high_list_primers <- tibble(
	primer_set = merged_high_list) %>% 
	mutate(
		Name_Forward = paste0(primer_set,"_f1"),
		Name_Reverse = paste0(primer_set,"_r1"))

write.csv(merged_high_list_primers, file = output_filename_remove, row.names = FALSE, quote = FALSE)


########
## 5- Diagnostic pdf
########

pdf(file=output_filename_diag_plot)

print(single_count_plot_pool)
print(boxplot_overall)
print(paired_line_plot)
print(heatmap_plot)
print(heatmap_plot_no_dil)
print(single_count_plot_pool_2)

dev.off()
