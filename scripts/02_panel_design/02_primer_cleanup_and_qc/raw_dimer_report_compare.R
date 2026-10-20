#module load R/4.4.2
packages <- c("dplyr", "data.table", "ggplot2","tidyverse","ggVennDiagram","gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
	install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)

# Distribution of Tm for MFE (how many happen before PCR temp and after, ones before will break apart)

# MFE delta G vs Pooler (if does not appear for one, set deltaG to 0):
# - shape = mfe, pooler, both
# - color = mfe_temp (blue = cold, red = high)

# boxplot deltaG mfe vs pooler (maybe separate, lower, higher, PCR temp maybe if points with boxplot diff color)
# boxplot only mfe and only pooler (add points, diff shape and colore for at PCR2 temp)

dimer_dir <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/"
raw_mfe_file <- paste(dimer_dir,"cisco_100wind_relax_primers_dimer_GTseq_param_mfe_parse.csv",sep="")
raw_pooler_file <- paste(dimer_dir,"pooler_bonds.csv",sep="")

raw_mfe <- fread(raw_mfe_file, header = TRUE)
raw_pooler <- fread(raw_pooler_file, header = TRUE)

# equivalent table
mfe_full_file <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/marker_real_primer_set.txt"
pooler_full_file <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/primers_primerpooler.tab"

col_names = c("ampliconID", "Forward", "Reverse")
mfe_full <- fread(mfe_full_file, header = FALSE, col.names = col_names)
pooler_full <- fread(pooler_full_file, header = FALSE, col.names = col_names)

PCR2_temp <- 57

## check size of each original table and resulting equivalent table
equiv_tb <- inner_join(
	mfe_full,pooler_full, by = "ampliconID",
	suffix = c(".mfe", ".pooler")) #change column names

### below need to make fct with compare_mfe_pooler
# 1st: make equivalent table in long format col mfe_primer, pooler_primer, F/Rev

long_equiv <- equiv_tb %>%
	pivot_longer(
		cols = c(Forward.mfe, Reverse.mfe, Forward.pooler, Reverse.pooler),
		names_to = c("Direction", ".value"),
		names_pattern = "(Forward|Reverse)\\.(mfe|pooler)"
		)

# First join: match Primer1 to pooler primers
pooler_annotated <- raw_pooler %>%
	left_join(long_equiv,
		by = c("primer1" = "pooler")) %>%
	rename(Primer1.mfe=mfe) %>%
	select(-c(Direction, ampliconID))

# Second join: match Primer2 to pooler primers
pooler_annotated <- pooler_annotated %>%
	left_join(long_equiv,
		by = c("primer2" = "pooler")) %>%
	rename(Primer2.mfe=mfe) %>%
	select(-c(Direction, ampliconID))

# 1st plot = histo temp (mfe only)
histo_temp <- ggplot(raw_mfe, aes(x = Tm)) +
			geom_histogram(aes(y = after_stat(density)), binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
			geom_density(color = "red", size = 1) +
			geom_vline(aes(xintercept = PCR2_temp), color = "red", linetype = "dashed", size = 1) +
			annotate("text", x = 54, y = 0.07, label = "PCR2", color = "red") +
			labs(title = "Distribution of MFE dimer Tm with Density Curve (all MFE dimers)",
				x = "Tm Range",
				y = "Density") +
			theme_bw() 

# 2nd plot = mfe vs pooler

raw_mfe <- raw_mfe %>% 
	mutate(
		pairs_name = paste(primer1, primer2, sep = "_")
		)


pooler_annotated <- pooler_annotated %>% 
	mutate(
		pairs_name = case_when( 
			paste(Primer1.mfe, Primer2.mfe, sep = "_") %in% raw_mfe$pairs_name ~ paste(Primer1.mfe, Primer2.mfe, sep = "_"),
			paste(Primer2.mfe, Primer1.mfe, sep = "_") %in% raw_mfe$pairs_name ~ paste(Primer2.mfe, Primer1.mfe, sep = "_"),
			TRUE ~ "pooler_specific" # Default 
		)
		) %>%
	select(-c(score,Tm))

#join into a new table keeping mfe and pooler specific bad pairs
comparison_tb <- full_join(raw_mfe,pooler_annotated, by = "pairs_name",
	suffix = c("_mfe", "_pooler")) %>%
	mutate(
		pairs_name = case_when( 
			is.na(deltaG_pooler) ~ "mfe_specific",
			TRUE ~ pairs_name
			)
		)%>%
	#select(pairs_name,primer1_mfe, primer2_mfe, deltaG_mfe, deltaG_pooler, Tm) 
	select(pairs_name,primer1_mfe, primer2_mfe, Primer1.mfe, Primer2.mfe, deltaG_mfe, deltaG_pooler, Tm) %>% rename(Primer1.mfe_pooler=Primer1.mfe, Primer2.mfe_pooler = Primer2.mfe)

#Warning message:
# In full_join(raw_mfe, pooler_annotated, by = "pairs_name", suffix = c("_mfe",  :
#   Detected an unexpected many-to-many relationship between `x` and `y`.
# ℹ Row 268 of `x` matches multiple rows in `y`.
# ℹ Row 292 of `y` matches multiple rows in `x`.
# ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.

# Convert NA (software specific dimers) into 0 (deltaG = 0 means no dimer) for plotting
# Add a type column to highlight the cases (software specific or shared)
comparison_tb <- comparison_tb %>%
	replace_na(list(deltaG_mfe = 0, deltaG_pooler = 0)) %>%
	mutate(
		type = case_when( 
			pairs_name == "mfe_specific" ~ "mfe_specific",
			pairs_name == "pooler_specific" ~ "pooler_specific",
			TRUE ~ "shared"
			))

dotplot <- ggplot(
	data = comparison_tb, aes(x = deltaG_mfe, y = deltaG_pooler, shape = type)
	) +
	# Map 'shape' to the factor type
	# Map 'fill' to the continuous variable temp for the gradient
	geom_point(
		aes(shape = type, fill = Tm),
		size = 4, # Adjust point size
		color = "black", # Set outline color
		stroke = 1) + # Adjust outline thickness
		#) +
	# Use scale_shape_manual to assign shapes 21-25 which have a fill aesthetic
  scale_shape_manual(values = c(21, 22, 23)) + # Use specific shapes that can be filled
  # Use a gradient scale function for the 'fill' aesthetic
  # scale_fill_gradient for a two-color gradient
  # scale_fill_gradientn for a multi-color gradient
  scale_fill_gradientn(
  	colors = c("blue", "yellow", "red"),
  	#values = c(0, 0.5, 1),              # Map these colors to specific data values (0, 0.5, 1)
  	values = scales::rescale(c(min(comparison_tb$Tm, na.rm = TRUE), 57, max(comparison_tb$Tm, na.rm = TRUE))), # Map colors to specific data values
  	breaks = c(min(comparison_tb$Tm), 57, max(comparison_tb$Tm)),              # Set breaks for the legend
  	labels = c("Low", "PCR2_57", "High"),
  	name = "mfe Tm") + # Customize the legend title
  #scale_fill_gradient(low = "blue", high = "red", name = "mfe Tm") +
	geom_abline(intercept = 0, slope = 1) +
	theme_bw() +
	ggtitle("DeltaG values pooler VS mfe")

# Boxplot focusing on only mfe and only pooler 
# (add points, diff shape and color for at PCR2 temp)
comparison_tb_spec <- comparison_tb %>%
	filter(pairs_name == "pooler_specific" | pairs_name == "mfe_specific") %>%
	pivot_longer(
		cols = c(deltaG_mfe, deltaG_pooler),
		#names_to = c("software", ".value"),
		#names_pattern = "(_mfe|_pooler)"
		names_to = "software",
    values_to = "deltaG"
		) %>%
	mutate(
		Tm_case = case_when( 
			Tm < 57 ~ "low",
			Tm > 57 ~ "high",
			Tm == 57 ~ "PCR2",
			TRUE ~ "NA"
		)) %>%
	# convert 0 to NA so will not be plotted
	mutate(deltaG = na_if(deltaG, 0))

boxplot_spec <- ggplot(
											comparison_tb_spec,
											aes(x = software, y = deltaG, color = software)
											) +
	geom_boxplot() +
	theme_bw() +
	geom_jitter(aes(shape = Tm_case),width = 0.2, alpha = 0.7) +
	labs(
		title = "DeltaG for program specific dimers",
		x = "Program",
		y = "DeltaG") +
	scale_shape_manual(values = c(15, 16, 17)) # Manually specify shapes

# Venn Diagram of primer pairs
#raw_mfe,pooler_annotated
# need a new pooler otherwise reduce to "pooler_specific"

pooler_annotated_venn <- raw_pooler %>%
	left_join(long_equiv,
		by = c("primer1" = "pooler")) %>%
	rename(Primer1.mfe=mfe) %>%
	select(-c(Direction, ampliconID))

# Second join: match Primer2 to pooler primers
pooler_annotated_venn <- pooler_annotated_venn %>%
	left_join(long_equiv,
		by = c("primer2" = "pooler")) %>%
	rename(Primer2.mfe=mfe) %>%
	select(-c(Direction, ampliconID))

pooler_annotated_venn <- pooler_annotated_venn %>% 
	mutate(
		pairs_name = case_when( 
			paste(Primer1.mfe, Primer2.mfe, sep = "_") %in% raw_mfe$pairs_name ~ paste(Primer1.mfe, Primer2.mfe, sep = "_"),
			paste(Primer2.mfe, Primer1.mfe, sep = "_") %in% raw_mfe$pairs_name ~ paste(Primer2.mfe, Primer1.mfe, sep = "_"),
			TRUE ~ paste(Primer1.mfe, Primer2.mfe, sep = "_") # Default 
		)
		) %>%
	select(-c(score,Tm))

venn_dt <- list(
	mfe = unique(raw_mfe$pairs_name),
	#pooler = pooler_annotated$pairs_name)
	pooler = unique(pooler_annotated_venn$pairs_name))

venn <- ggVennDiagram(venn_dt) + 
	ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2)) +
	ggtitle("Number of Dimers identified")

# temp distribution for the ones shared between programs
histo_temp_shared <- ggplot(
					(comparison_tb %>%
					filter(pairs_name != "pooler_specific" & pairs_name != "mfe_specific") ), 
				aes(x = Tm)) +
			geom_histogram(aes(y = after_stat(density)), binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
			geom_density(color = "red", size = 1) +
			geom_vline(aes(xintercept = PCR2_temp), color = "red", linetype = "dashed", size = 1) +
			annotate("text", x = 54, y = 0.07, label = "PCR2", color = "red") +
			labs(title = "Distribution of MFE dimer Tm with Density Curve (shared dimers)",
				x = "Tm Range",
				y = "Density") +
				ylim(0, 0.08) + # Sets the y-axis from 0 to 0.08 so comparable with 1st graph
			theme_bw() 

# to add:
# mfe rest deltaG vs temp, shared with pooler in diff shape
delta_temp_plot <- ggplot(
	data = (comparison_tb %>%
					filter(pairs_name != "pooler_specific")), 
	aes(x = Tm, y = deltaG_mfe, shape = type)
	) +
	geom_point(
		aes(shape = type, fill = type),
		size = 4, # Adjust point size
		color = "black", # Set outline color
		stroke = 1) +
	scale_shape_manual(values = c(21, 22)) + 
	geom_vline(aes(xintercept = PCR2_temp), color = "red", linetype = "dashed", size = 1) +
	labs(title = "MFE dimer deltaG vs Tm") +
	theme_bw()



# Save pdf as an overlow diagnostic (not each plot separately)

# Open the PDF device
pdf("/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/compare_raw_mfe_pooler_diagnostic_plots.pdf", width = 8, height = 10)

# Arrange and print the plots
#grid.arrange(venn,dotplot, heatmap, nrow = 3) # Arrange in 2 columns
# Move to a new page
#plot.new()
#print(histo_temp)
#print(histo_temp_shared)
grid.arrange(histo_temp,histo_temp_shared, nrow = 1) # Arrange in 2 columns
print(dotplot)
print(boxplot_spec)
print(venn)
print(delta_temp_plot)
# Close the PDF device
dev.off()

#######################
#### PCR test
#######################

# Output: .xls table with following columns
## - PCR # (keep track of number of reaction)
## - Primer 1 name
## - Primer 2 name
## - DeltaG MFEprimer
## - Temp MFEprimer
## - DeltaG primerpooler
## - Dimer expectations: 
### * high (if temp == PCR2 & DeltaG MFE =< -9 & primerpooler =< -3)
### * Medium (if temp / PCR2 & DeltaG MFE =< -9 & primerpooler =< -3)
### * Medium (if temp == PCR2 & -9 < DeltaG MFE <0 & -3 < primerpooler <0)
### * Low (if temp == PCR2 & DeltaG MFE =<-9 & primerpooler == 0)
### * Low (DeltaG MFE ==0 & primerpooler =< -3)
### * Low (if temp == PCR2 & -9 < DeltaG MFE <0 & primerpooler == 0)
### * Low (DeltaG MFE ==0 & -3 < primerpooler < 0)
### * None (DeltaG MFE ==0, not detected by both programs - select 2 bad primers that should not interact to save good primers)
## - Dimer PCR (empty, to be filled)

PCR_test_tb <- comparison_tb %>%
	#select(-pairs_name) %>%
	mutate(
		Dimer_expectations = case_when( 
			type == "shared" & round(Tm) == PCR2_temp & deltaG_mfe <= -9 & deltaG_pooler <= -3 ~ "High",
			type == "shared" & round(Tm) != PCR2_temp & deltaG_mfe <= -9 & deltaG_pooler <= -3 ~ "Medium",
			type == "shared" & ( deltaG_mfe > -9 | deltaG_pooler > -3 ) ~ "Medium",
			type == "mfe_specific" | type == "pooler_specific" ~ "Low",
			TRUE ~ "NA" # Default 
			#TRUE ~ "NA"
		),
		limitations = case_when(
			Dimer_expectations == "High" ~ "none",
			Dimer_expectations == "Medium" & round(Tm) != PCR2_temp ~ "Tm",
			#Dimer_expectations == "Medium" & Tm == PCR2_temp ~ "Tm",
			Dimer_expectations == "Medium" & round(Tm) == PCR2_temp & (deltaG_mfe > -9 | deltaG_pooler > -3)~ "highDeltaG",
			type == "mfe_specific" & round(Tm) == PCR2_temp ~ "mfe_specific",
			type == "mfe_specific" & round(Tm) != PCR2_temp ~ "mfe_specific_Tm",
			type == "pooler_specific" ~ "pooler_specific",
			TRUE ~ "other" # Default 
			)#,
		#Dimer_expectations = as.factor(Dimer_expectations),
		#limitations = as.factor(limitations)
		) %>%
	mutate(
		Primer_1 = case_when( 
			is.na(primer1_mfe) ~ Primer1.mfe_pooler,
			TRUE ~ primer1_mfe), # Default,
			Primer_2 = case_when( 
			is.na(primer2_mfe) ~ Primer2.mfe_pooler,
			TRUE ~ primer2_mfe # Default
		)) %>%
	dplyr::group_by(Dimer_expectations,limitations) %>%
	dplyr::slice_head(n=5) %>%
	#dplyr::arrange(deltaG_mfe)%>%
	#dplyr::slice(c(1, n())) %>%
	#dplyr:: filter(row_number() == 1 | row_number() == n())
	#dplyr::bind_rows(dplyr::slice_max(., order_by = deltaG_mfe, n = 5), dplyr::slice_min(., order_by = deltaG_mfe, n = 5)) %>%
	ungroup() %>% select(Primer_1,Primer_2, deltaG_mfe:limitations)

PCR_test_tb <-	PCR_test_tb %>% mutate(rownames=rownames(PCR_test_tb), PCR_reaction=paste("PCR_",rownames,sep=""),
	Dimer_PCR="")

print(PCR_test_tb)
head(PCR_test_tb)

# need to get "bad" primers that interract ok together only
bad <- PCR_test_tb %>%
  drop_na(Primer_1) %>% 
  #head() %>% 
  select(Primer_1) %>%
  unique()%>%
  dplyr::slice_head(n=5)
  

good_interaction <- PCR_test_tb %>% filter( ! Primer_1 %in% bad$Primer_1 & ! Primer_2 %in% bad$Primer_1) %>%
	drop_na(Primer_1,Primer_2) 

#now we can associate bad with either primer1_mfe or primer2_mfe 
primer_friends <- sample(unique(c(good_interaction$Primer_1,good_interaction$Primer_2)),5, replace = FALSE)

# now we have our good sets of bad primers just need to add to table
# create a temp table with all the colnames then mutate the specific columns with the vectors
# then rbind

empty_tb <- data.frame(
  matrix(nrow = length(bad$Primer_1), ncol = length(colnames(PCR_test_tb)))
) 

colnames(empty_tb) <- colnames(PCR_test_tb)
empty_tb <- empty_tb %>% mutate(
	Primer_1 = bad$Primer_1,
  Primer_2 = primer_friends,
  Dimer_expectations = "none")

str(empty_tb)
PCR_test_tb <- rbind(PCR_test_tb,empty_tb)

write.table(
	PCR_test_tb,
	"/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/PCR_test_dimers.csv", 
	sep = ",",row.names=FALSE, quote = FALSE
	)

# Get the selected primers in good ordering format
primer <- fread("/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/tm57/cisco_100wind_relax_primers_ordering_format_all_dimers.csv",sep=",")

primer_test <- primer %>% filter(Name %in% c(PCR_test_tb$Primer_1,PCR_test_tb$Primer_2) )

#check
PCR_test_tb %>% filter( !PCR_test_tb$Primer_1 %in% primer_test$Name) 
PCR_test_tb %>% filter( !PCR_test_tb$Primer_2 %in% primer_test$Name) 

primer_test <- primer_test %>% select(Name,Sequence)

write.table(
	primer_test,
	"/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/cisco_100wind_test_dimers.csv", 
	sep = ",",row.names=FALSE, quote = FALSE
	)