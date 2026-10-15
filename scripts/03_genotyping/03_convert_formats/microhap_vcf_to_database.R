###############################################################################
# microhap_vcf_to_database.R
#
# Purpose:
# Convert a microhap VCF (GTseq projected SNPs or GATK calls) into a
# database-ready SNP matrix.
#
# Input:
#   VCF file given with --vcf (run with --help for all options)
#
# Output (file_suffix is appended to each name):
#   microhap_database.csv
#   microhap_database_bloater.csv
#   locus_summary.tsv
#   filtering_summary.tsv
#   microhap_pca.png
#   variants_per_locus_histogram.png
#   sample_GTmissingness.png
#
# Output format:
#
# GEN_ID      57853_53_1   57853_53_2   12345_17_1
# BLNE-18-001 CC           GG           TT
# BLNE-23-003 CC           GG           CT
# 604-0005    CG           GT           TT
#
# Notes:
#   - Each VCF record is already a SNP.
#   - SNPs are grouped by LOCUS. How LOCUS is defined depends on vcf_source:
#
#       "gtseq" : parsed from the ID field
#                 57853_53:SNP:40368066 -> 57853_53
#
#       "gatk"  : taken from the CHROM column (ID is usually "." in GATK VCFs)
#                 e.g. CHROM 57853_53 -> 57853_53
#
#   - SNPs are numbered by genomic position (POS) within each locus.
#   - Missing calls are recorded as "0".
#   - Sample names are cleaned:
#
#       001_BLNE-18-001_10489050_GY8VH
#       -> BLNE-18-001
#
#       097_CT367_604-0005_10489050_GY8VH
#       -> 604-0005
#
#       096_NEFC_BloaterPlate_NTC_10489050_GY8VH
#       -> NEFC_BloaterPlate_NTC
#
###############################################################################
#module load R/4.4.2
library(vcfR)
library(dplyr)
library(tibble)
library(ggplot2)

###############################################################################
# Command-line arguments
###############################################################################

usage <- function(status = 0){
	cat(
"
Usage:
  Rscript microhap_vcf_to_database.R --vcf <file.vcf> [options]

Required:
  --vcf <path>            Input VCF file (GTseq microhap or GATK)

Optional:
  --out_dir <dir>         Output directory (default: current directory)
  --vcf_source <type>     'gtseq' or 'gatk' (default: gtseq)
                            gtseq: locus parsed from the ID field
                            gatk : locus taken from the CHROM column
  --ind_missing <0-1>     Max missingness per individual (default: 0.80)
  --site_missing <0-1>    Max missingness per site (default: 0.80)
  --suffix <string>       Suffix added to output file names (default: _relax)
  -h, --help              Show this message

Arguments can be given as '--key value' or '--key=value'.

Example:
	module load R/4.4.2
  Rscript microhap_vcf_to_database.R \\
    --vcf /path/to/microhap_variants.vcf \\
    --out_dir /path/to/output/ \\
    --vcf_source gatk --suffix _gatk
\n", sep = "")
	quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if(any(args %in% c("-h", "--help"))) usage(0)

# Defaults (vcf has none and is required)
opts <- list(
	vcf          = NA_character_,
	out_dir      = ".",
	vcf_source   = "gtseq",
	ind_missing  = "0.80",
	site_missing = "0.80",
	suffix       = "_relax"
)

i <- 1
while(i <= length(args)){

	a <- args[i]

	if(!grepl("^--", a)){
		cat("Unexpected argument: ", a, "\n", sep = "")
		usage(1)
	}

	if(grepl("=", a, fixed = TRUE)){
		key <- sub("^--([^=]+)=.*$", "\\1", a)
		val <- sub("^--[^=]+=", "", a)
		i <- i + 1
	} else {
		key <- sub("^--", "", a)
		if(i == length(args)){
			cat("Missing value for --", key, "\n", sep = "")
			usage(1)
		}
		val <- args[i + 1]
		i <- i + 2
	}

	key <- gsub("-", "_", key)

	if(!key %in% names(opts)){
		cat("Unknown argument: --", key, "\n", sep = "")
		usage(1)
	}

	opts[[key]] <- val
}

if(is.na(opts$vcf)){
	cat("Error: --vcf is required\n")
	usage(1)
}

vcf_path               <- opts$vcf
out_dir                <- opts$out_dir
vcf_source             <- opts$vcf_source
ind_missing_threshold  <- suppressWarnings(as.numeric(opts$ind_missing))
site_missing_threshold <- suppressWarnings(as.numeric(opts$site_missing))
file_suffix            <- opts$suffix

if(!file.exists(vcf_path))
	stop("VCF file not found: ", vcf_path)

if(!vcf_source %in% c("gtseq", "gatk"))
	stop("--vcf_source must be 'gtseq' or 'gatk', not '", vcf_source, "'")

if(is.na(ind_missing_threshold) ||
   ind_missing_threshold < 0 || ind_missing_threshold > 1)
	stop("--ind_missing must be a number between 0 and 1")

if(is.na(site_missing_threshold) ||
   site_missing_threshold < 0 || site_missing_threshold > 1)
	stop("--site_missing must be a number between 0 and 1")

# Output paths are built with paste0(out_dir, ...), so ensure trailing slash
if(!grepl("/$", out_dir)) out_dir <- paste0(out_dir, "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Read VCF
###############################################################################

vcf <- read.vcfR(
	vcf_path,
	verbose = FALSE
)

###############################################################################
# Extract fixed fields and genotypes
###############################################################################

fix <- as.data.frame(
	getFIX(vcf),
	stringsAsFactors = FALSE
)

# getFIX returns everything as character; POS must be numeric so that
# SNPs are ordered by position (otherwise "100" sorts before "20")
fix$POS <- as.numeric(fix$POS)

gt <- extract.gt(
	vcf,
	element = "GT"
)

###############################################################################
# Extract locus name
###############################################################################

if(vcf_source == "gtseq"){

	# Extract locus from ID
	#
	# 57853_53:SNP:40368066
	# -> 57853_53

	fix$LOCUS <- sub(
		":.*$",
		"",
		fix$ID
	)

} else if(vcf_source == "gatk"){

	# GATK VCFs have no locus information in ID (usually "."),
	# so the locus is the CHROM (reference sequence) name.
	# POS is used below to order and number SNPs within each locus.

	fix$LOCUS <- fix$CHROM
}

###############################################################################
# Number SNPs within loci
###############################################################################

# Keep track of the original row order so gt can be re-ordered to match fix
# after sorting (fix is sorted by LOCUS then POS, gt is still in VCF order)
fix$ROW_IDX <- seq_len(nrow(fix))

fix <- fix %>%
	group_by(LOCUS) %>%
	arrange(POS, .by_group = TRUE) %>%
	mutate(
		SNP_NUM = row_number(),
		SNP_ID = paste0(LOCUS, "_", SNP_NUM)
	) %>%
	ungroup()

gt <- gt[fix$ROW_IDX, , drop = FALSE]

###############################################################################
# Clean sample names
###############################################################################

clean_sample_id <- function(x){

	# Remove leading well number
	x <- sub(
		"^[0-9]+_",
		"",
		x
	)

	# Remove run ID and flowcell suffix
	x <- sub(
		"_[0-9]+_[^_]+$",
		"",
		x
	)

	# Remove Cisco panel prefix
	x <- sub(
		"^CT367_",
		"",
		x
	)

	x
}

sample_ids <- clean_sample_id(
	colnames(gt)
)

###############################################################################
# Missingness filtering
###############################################################################

#is_missing <- gt == "./." | gt == ".|."

#--------------------------------------------------
# Remove individuals
#--------------------------------------------------

#ind_missing <- colMeans(is_missing)
n_ind_start <- ncol(gt) # saving for summary output
ind_missing <- colMeans(is.na(gt)) #missing = NA with vcfR

keep_ind <- ind_missing <= ind_missing_threshold

n_ind_removed_missing <- sum(!keep_ind)

gt <- gt[, keep_ind, drop = FALSE]

sample_ids <- sample_ids[keep_ind]

#--------------------------------------------------
# Remove variants AFTER individual filtering
#--------------------------------------------------

#site_missing <- rowMeans(is_missing)
n_site_before_missing <- nrow(gt)
site_missing <- rowMeans(is.na(gt))

keep_site <- site_missing <= site_missing_threshold
n_site_removed_missing <- sum(!keep_site)

gt <- gt[keep_site, , drop = FALSE]

fix <- fix[keep_site, ]

###############################################################################
# Remove invariant markers
###############################################################################

n_site_before_invariant <- nrow(gt)
n_genotypes <- apply(
	gt,
	1,
	function(x){

		#x <- x[!(x %in% c("./.", ".|."))]
		x <- x[!is.na(x)]

		# Normalise phasing and allele order so 0/1, 1/0, 0|1 and 1|0 count as
		# the same genotype (otherwise a site with only heterozygotes written
		# in different ways is wrongly kept as "variable")
		x <- gsub("|", "/", x, fixed = TRUE)
		x <- vapply(
			strsplit(x, "/", fixed = TRUE),
			function(a) paste(sort(a), collapse = "/"),
			character(1)
		)

		length(unique(x))
	}
)

keep_site <- n_genotypes > 1
n_site_removed_invariant <- sum(!keep_site)

gt <- gt[keep_site, ]
fix <- fix[keep_site, ]

###############################################################################
# Remove multiallelic variants
###############################################################################

keep_site <- !grepl(",", fix$ALT)

n_site_removed_multiallelic <- sum(!keep_site)

gt  <- gt[keep_site, , drop = FALSE]
fix <- fix[keep_site, ]

###############################################################################
# PCA
###############################################################################

# Convert genotypes to allele counts
# 0/0 = 0
# 0/1 = 1
# 1/1 = 2
# ./. = NA

gt_num <- matrix(
	NA_real_,
	nrow = ncol(gt),
	ncol = nrow(gt),
	dimnames = list(
		sample_ids,
		fix$SNP_ID
	)
)

for(i in seq_len(nrow(gt))) {

	g <- gt[i, ]

	gt_num[, i] <- ifelse(
		g %in% c("0/0", "0|0"), 0,
		ifelse(
			g %in% c("0/1", "1/0", "0|1", "1|0"), 1,
			ifelse(
				g %in% c("1/1", "1|1"), 2,
				NA
			)
		)
	)
}

# Impute missing values using marker means
for(i in seq_len(ncol(gt_num))) {

	m <- mean(gt_num[, i], na.rm = TRUE)

	gt_num[is.na(gt_num[, i]), i] <- m
}

###############################################################################
# Species assignment from sample names
###############################################################################

species <- ifelse(
	grepl("^BLNE-", sample_ids),
	"Bloater",
	ifelse(
		grepl("^NEFC", sample_ids),
		"Control",
		"Cisco"
	)
)

###############################################################################
# Remove invariant markers
###############################################################################

# vars <- apply(
# 	gt_num,
# 	2,
# 	var,
# 	na.rm = TRUE
# )

# keep <- vars > 0

# summary(vars)
# table(is.na(vars))

# gt_num_pca <- gt_num[, keep]

# sds <- apply(gt_num_pca, 2, sd) #remove if standard dev is 0 between 2 SNPS

# keep <- !is.na(sds)

# gt_num_pca <- gt_num_pca[, keep]

### sanity check on invariants

sds <- apply(gt_num, 2, sd)

bad <- which(is.na(sds) | sds == 0)

if(length(bad) > 0){

	cat("\nERROR: ", length(bad), " SNP(s) have zero variance in the PCA matrix.\n",
		sep = "")
	cat("Raw genotype counts for the first few:\n\n")

	for(b in head(bad, 10)){
		cat(colnames(gt_num)[b], "\n")
		print(table(gt[b, ], useNA = "ifany"))
		cat("\n")
	}

	stop("Zero-variance SNPs found; see genotype tables above")
}

###############################################################################
# Run PCA
###############################################################################

pca <- prcomp(
	gt_num,
	center = TRUE,
	scale. = TRUE
)

pca_df <- data.frame(
	GEN_ID = sample_ids,
	Species = species,
	PC1 = pca$x[,1],
	PC2 = pca$x[,2]
)

###############################################################################
# Plot PCA
###############################################################################

pca <- ggplot(
	pca_df,
	aes(
		PC1,
		PC2,
		colour = Species
	)
) +
	geom_point(size = 3, alpha = 0.8) +
	theme_bw() +
	labs(
		title = "PCA of Microhap Variants"
	)

ggsave(
	paste0(out_dir,"microhap", file_suffix, "_pca.png"),
	pca,
	width = 8,
	height = 6,
	dpi = 300
)

# ###############################################################################
# # Save sample manifest
# ###############################################################################

# sample_manifest <- data.frame(
# 	ORIGINAL_NAME = colnames(gt),
# 	GEN_ID = sample_ids,
# 	stringsAsFactors = FALSE
# )

# write.table(
# 	sample_manifest,
# 	"sample_manifest.tsv",
# 	sep = "\t",
# 	quote = FALSE,
# 	row.names = FALSE
# )

###############################################################################
# Convert GT calls to allele strings
###############################################################################

convert_gt <- function(gt_call, ref, alt){

	if(is.na(gt_call) || gt_call == "./." || gt_call == ".|.")
		return("0")

	if(gt_call %in% c("0/0", "0|0"))
		return(paste0(ref, ref))

	# if(gt_call %in% c(
	# 	"0/1", "1/0", ##### NEED TO SEPARATE
	# 	"0|1", "1|0"
	# ))
	# 	return(paste0(ref, alt))

	if(gt_call %in% c(
		"0/1", "0|1" ##### NEED TO SEPARATE
	))
		return(paste0(ref, alt))

	if(gt_call %in% c(
		"1/0", "1|0" ##### NEED TO SEPARATE
	))
		return(paste0(alt, ref))

	if(gt_call %in% c("1/1", "1|1"))
		return(paste0(alt, alt))

	return("0")
}

###############################################################################
# Create SNP matrix
###############################################################################

out <- matrix(
	NA_character_,
	nrow = ncol(gt),
	ncol = nrow(gt),
	dimnames = list(
		sample_ids,
		fix$SNP_ID
	)
)

for(i in seq_len(nrow(gt))){

	out[, i] <- vapply(
		gt[i, ],
		convert_gt,
		character(1),
		ref = fix$REF[i],
		alt = fix$ALT[i]
	)
}

###############################################################################
# Convert to data frame
###############################################################################

out <- as.data.frame(
	out,
	stringsAsFactors = FALSE
)

out <- rownames_to_column(
	out,
	var = "GEN_ID"
)

###############################################################################
# Write database file
###############################################################################

# write.table(
# 	out,
# 	paste0(out_dir,"microhap_database.tsv"),
# 	sep = "\t",
# 	quote = FALSE,
# 	row.names = FALSE
# )

write.csv(
	out,
	paste0(out_dir, "microhap_database", file_suffix, ".csv"),
	row.names = FALSE,
	quote = FALSE
)

###############################################################################
# Create Bloater-only database export
###############################################################################

out_bloater <- subset(
	out,
	grepl("^BLNE-", GEN_ID)
)

# write.table(
# 	out_bloater,
# 	paste0(out_dir, "microhap_database_bloater.tsv"),
# 	sep = "\t",
# 	quote = FALSE,
# 	row.names = FALSE
# )

write.csv(
	out_bloater,
	paste0(out_dir, "microhap_database_bloater", file_suffix, ".csv"),
	row.names = FALSE,
	quote = FALSE
)

# Use the same species assignment as the PCA (rows of out follow sample_ids)
n_bloater <- sum(species == "Bloater")
n_cisco   <- sum(species == "Cisco")
n_control <- sum(species == "Control")

stopifnot(n_bloater + n_cisco + n_control == nrow(out))

###############################################################################
# SNPs per locus
###############################################################################

locus_summary <- fix %>%
	group_by(LOCUS) %>%
	summarise(
		N_VARIANTS = n(),
		.groups = "drop"
	)

write.table(
	locus_summary,
	paste0(out_dir, "locus_summary", file_suffix, ".tsv"),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE
)

locus_hist <- ggplot(
	locus_summary,
	aes(x = N_VARIANTS)
) +
	geom_histogram(
		binwidth = 1,
		color = "black",
		fill = "steelblue"
	) +
	theme_bw() +
	labs(
		title = "Variants per Microhap Locus",
		x = "Number of Variants",
		y = "Number of Loci"
	)

ggsave(
	paste0(out_dir, "variants_per_locus", file_suffix, "_histogram.png"),
	locus_hist,
	width = 8,
	height = 6,
	dpi = 300
)

n_loci <- n_distinct(fix$LOCUS)

mean_var_per_locus <- round(
	mean(locus_summary$N_VARIANTS),
	2
)

max_var_per_locus <- max(locus_summary$N_VARIANTS)

median_var_per_locus <- median(locus_summary$N_VARIANTS)

###############################################################################
# Filtering Summary
###############################################################################

cat("\n")
cat("Filtering summary\n")
cat("-----------------\n")
cat("Individuals removed (missingness): ",
	n_ind_removed_missing, "\n", sep = "")
cat("Sites removed (missingness): ",
	n_site_removed_missing, "\n", sep = "")
cat("Sites removed (invariant): ",
	n_site_removed_invariant, "\n", sep = "")
cat("Sites removed (multiallelic): ",
	n_site_removed_multiallelic, "\n", sep = "")
cat("Final individuals: ", ncol(gt), "\n", sep = "")
cat("Final Bloaters: ", n_bloater, "\n", sep = "")
cat("Final Ciscoes: ", n_cisco, "\n", sep = "")
cat("Final sites: ", nrow(gt), "\n", sep = "")
cat("Final loci: ", n_loci, "\n", sep = "")
cat("Max SNPs per locus: ", max_var_per_locus, "\n", sep = "")
cat("Mean SNPs per locus: ", mean_var_per_locus, "\n", sep = "")
cat("Median SNPs per locus: ", median_var_per_locus, "\n", sep = "")
cat("\n")

# filtering_summary <- data.frame(
#   Metric = c(
# 	"Individuals_initial",
# 	"Individuals_removed_missingness",
# 	"Individuals_final",
# 	"Sites_initial",
# 	"Sites_removed_missingness",
# 	"Sites_removed_invariant",
# 	"Sites_final"
#   ),
#   Value = c(
# 	n_ind_start,
# 	n_ind_removed_missing,
# 	ncol(gt),
# 	n_site_start,
# 	n_site_removed_missing,
# 	n_site_removed_invariant,
# 	nrow(gt)
#   )
# )

filtering_summary <- data.frame(
	Metric = c(
		"Individuals_initial",
		"Individuals_removed_missingness",
		"Individuals_final",
		"Bloaters_final",
		"Cisco_final",
		"Controls_final",
		"Sites_initial",
		"Sites_removed_missingness",
		"Sites_removed_invariant",
		"Sites_removed_multiallelic",
		"Sites_final",
		"Loci_final",
		"Max_sites_per_locus",
		"Mean_sites_per_locus",
		"Median_sites_per_locus"
	),
	Value = c(
		n_ind_start,
		n_ind_removed_missing,
		ncol(gt),
		n_bloater,
		n_cisco,
		n_control,
		n_site_before_missing,
		n_site_removed_missing,
		n_site_removed_invariant,
		n_site_removed_multiallelic,
		nrow(gt),
		n_loci,
		max_var_per_locus,
		mean_var_per_locus,
		median_var_per_locus
	)
)

write.table(
	filtering_summary,
	paste0(out_dir, "filtering", file_suffix, "_summary.tsv"),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE
)

###############################################################################
# Sample-level QC
###############################################################################

missing_by_sample <- data.frame(
	GEN_ID = out$GEN_ID,
	MissingPct = rowMeans(out[, -1] == "0") * 100
)

# write.table(
# 	missing_by_sample,
# 	"sample_missingness.tsv",
# 	sep = "\t",
# 	quote = FALSE,
# 	row.names = FALSE
# )

###############################################################################
# Plot sample missingness
###############################################################################

GTmissingness_plot <- ggplot(
	missing_by_sample,
	aes(
		x = reorder(GEN_ID, MissingPct),
		y = MissingPct
	)
) +
	geom_col() +
	coord_flip() +
	labs(
		title = "Missing Data by Sample",
		x = "Sample",
		y = "Percent Missing"
	) +
	theme_bw() +
	theme(
		axis.text.y = element_text(size = 4)
		)

ggsave(
	paste0(out_dir,"sample_GTmissingness", file_suffix, ".png"),
	GTmissingness_plot,
	width = 8,
	height = 10,
	dpi = 300
)

###############################################################################
# Summary
###############################################################################

cat("\n")
cat("Samples:", nrow(out), "\n")
cat("SNPs:", ncol(out) - 1, "\n")
cat("Database file written: ",
	paste0(out_dir, "microhap_database", file_suffix, ".csv"), "\n", sep = "")
cat("Bloater database file written: ",
	paste0(out_dir, "microhap_database_bloater", file_suffix, ".csv"), "\n", sep = "")
cat("\n")
