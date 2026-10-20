#module load R/4.4.2
packages <- c("dplyr", "data.table", "ggplot2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)

# Output directorys
out_eval <- "/local/workdir/cauret/gtseq_cisco/seq_panel_1/eval/"

# List all read count files
count_files <- list.files(path = out_eval, pattern = "\\.txt$", full.names = TRUE)

# Read each CSV file into a list of data frames
list_of_dfs <- lapply(count_files, read.csv)

# Combine all data frames into a single data frame
count_table <- do.call(rbind, list_of_dfs) %>%
  select(-file)

# View the first few rows of the combined table
head(count_table)

# remove markers where 0 read for all individuals
# check how many markers left for plotting purposes


filtered_table <- count_table %>%
  group_by(primer_set) %>%
  summarise(total_primers_sum = sum(total_primers)) %>%
  filter(total_primers_sum > 0) %>%
  inner_join(count_table, by = "primer_set")

length(unique(filtered_table$primer_set))
length(unique(count_table$primer_set))

#####
# 1- Variation of primer counts between individuals for each set
#####

# Since. >600 primer_sets, make groups

filtered_table <- filtered_table %>%
  mutate(primer_bin = case_when(
    total_primers >= 0 & total_primers < 5000 ~ "Low",
    total_primers >= 5000 & total_primers < 10000 ~ "Medium",
    total_primers >= 10000 ~ "High"
  ), primer_bin = factor(primer_bin, levels = c("Low", "Medium", "High"))
)

# Check how many data points are in each bin
table(filtered_table$primer_bin)

# Create boxplot (not enough variance for violin)
box_count <- ggplot(filtered_table, aes(x = primer_bin, y = total_primers)) +
  #geom_violin(trim = FALSE, fill = "darkblue", color = "black", na.rm = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_violin(fill = "darkblue", na.rm = TRUE) +
  geom_jitter(width = 0.01, alpha = 0.05, color = "lightblue") +
  annotate("text", x = 1, y = max(count_table$total_primers) * 1.05, label = "0–4999", size = 4) +
  annotate("text", x = 2, y = max(count_table$total_primers) * 1.05, label = "5000–9999", size = 4) +
  annotate("text", x = 3, y = max(count_table$total_primers) * 1.05, label = "10000+", size = 4) +
  labs(title = "Primer Set Distribution by Bin",
       x = "Primer Bin",
       y = "Total Primers") +
  theme_minimal()
ggsave(box_count, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/box_bin_read_count.png")


# Create violin plot with boxplot overlay
violin_count <- ggplot(filtered_table, aes(x = primer_set, y = total_primers)) +
  #geom_violin(trim = FALSE, fill = "lightblue", color = "black") +
  geom_violin(fill = "lightblue", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ primer_bin, scales = "free_x") +
  #geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(title = "Variation of Total Primers Across Primer Sets",
       x = "Primer Set",
       y = "Total Primers") +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x = element_blank()
    )

ggsave(violin_count, 
  width = 20,
  height = 6,
  units = "in", 
  filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/violin_read_count.png")

# Top N Primer Sets by Total Primers
top_sets <- filtered_table %>%
  group_by(primer_set) %>%
  summarise(total_primers_sum = sum(total_primers)) %>%
  arrange(desc(total_primers_sum)) %>%
  slice_head(n = 30) %>%
  pull(primer_set)

top_sets_table <- filtered_table %>%
  filter(primer_set %in% top_sets)

violin_count_top_sets <- ggplot(top_sets_table, aes(x = primer_set, y = total_primers)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.05) +
  labs(title = "Top 30 Primer Sets by Total Primers",
       x = "Primer Set", y = "Total Primers") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(violin_count_top_sets, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/violin_read_count_top_sets.png")


## Focus on reads with 2 primers
# expected: F/R_rc F_rc/R
# not_expected: F/R F_rc/R_rc

#remove if no reads for those category in all individuals

filtered_table <- filtered_table %>%
  mutate(
    total_single = only_forward + only_forward_rc + only_reverse + only_reverse_rc,
    total_expected_pair = forward_reverse_rc + forward_rc_reverse,
    total_wrong_orientation_pair = forward_reverse + forward_rc_reverse_rc) #wrong orientation

#mean by categorie for each primer_set, then count number of primer set for each cat
filtered_table_sum_cat <- filtered_table %>% 
  group_by(primer_set) %>%
  summarise(mean_single = mean(total_single),
    mean_expected_pair = mean(total_expected_pair),
    mean_wrong_orientation_pair = mean(total_wrong_orientation_pair))

filtered_table_sum_cat %>%
  filter(mean_expected_pair > 0) %>%
  #filter(mean_expected_pair > 5) %>% #tally - 0
  tally() #118

filtered_table_sum_cat %>%
  filter(mean_expected_pair ==0) %>%
  tally() #508

# heatmap with markers at bottom and ind on side, color = numb reads
#maybe add colm that has any combination of forw and rev?
#need to remove file ext from ind, change orientation of name on x axis 
hp_read_count <- ggplot(results_df, 
  aes(x = primer_set, y = individuals, fill = total_primers,color = "black")) +
  geom_tile(color = "black",
            #color = "white",
            #lwd = 1.5,
            linetype = 1
            ) +
  scale_fill_gradient(low = "white", high = "red") +
  #coord_fixed() +
  theme(axis.text.x = element_blank()#, # Remove x-axis tick labels
        #axis.text.y = element_blank()
        )#,
        #axis.title.x = element_blank(), # Remove x-axis title
        #axis.title.y = element_blank()) # Remove y-axis title # Remove y-axis tick labels
#ggsave(hp_read_count, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/hp_read_count.png")
ggsave(hp_read_count,
  width = 20,
  height = 6,
  units = "in", 
  filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/hp_read_count.png")


## check if some individuals worked
min_reads = 2
#min_reads = 5

filtered_table %>% filter(forward_reverse > min_reads | 
                          forward_rc_reverse > min_reads | 
                          forward_reverse_rc > min_reads | 
                          forward_rc_reverse_rc > min_reads ) %>%
                  #select(forward_reverse:individuals)
                  #select(primer_set,forward_reverse_rc:individuals)
                  select(primer_set,individuals)

# # A tibble: 6 × 5
#   primer_set forward_reverse_rc forward_rc_reverse_rc total_primers individuals 
#   <chr>                   <int>                 <int>         <int> <chr>       
# 1 2329125_58                  0                     6         13785 CT_cART_465…
# 2 2444222_38                  7                     0             7 CT_cHOY_122…
# 3 2444222_38                 11                     0            11 CT_cHOY_123…
# 4 793800_58                   0                     6         13374 CT_Tow1_noS…
# 5 793800_58                   0                     6         12012 CT_cART_465…
# 6 793800_58                   0                     6         12458 CT_cART_13_…

# # So the only "success" would be if min sert to 5 for expected good combination (forward_rc_reverse, forward_reverse_rc)
# 2 2444222_38 CT_cHOY_122_10483695_M53J6_R1         
# 3 2444222_38 CT_cHOY_123_10483695_M53J6_R1         


OK_sub <- filtered_table %>% filter(#forward_reverse > min_reads | 
                          forward_rc_reverse >= min_reads | 
                          forward_reverse_rc >= min_reads #| 
                          #forward_rc_reverse_rc > min_reads 
                          ) #%>%
                  #select(forward_reverse:individuals)
                  #select(primer_set,forward_reverse_rc:individuals)
                  #select(primer_set,individuals)
OK_sub %>% select(primer_set,individuals, forward_reverse_rc, forward_rc_reverse)

length(unique(OK_sub$primer_set))

# Group by 'cyl' and 'gear', then summarise with a count
result_summarise <- OK_sub %>%
  group_by(individuals) %>%
  summarise(number_of_lines = n(), .groups = 'drop') # .groups = 'drop' removes grouping after the operation

# View the result
print(result_summarise)
result_summarise %>% filter(number_of_lines >1)

# result_summarise %>% filter(number_of_lines >1)
# # A tibble: 12 × 2
#    individuals                            number_of_lines
#    <chr>                                            <int>
#  1 CT_Tow1_noSQP_300_9_10483695_M53J6_R1                2
#  2 CT_Tow1_noSQP_455_21_10483695_M53J6_R1               2
#  3 CT_cART_13_10483695_M53J6_R1                         4
#  4 CT_cART_466_23_10483695_M53J6_R1                     2
#  5 CT_cART_F2Q_10483695_M53J6_R1                        3
#  6 CT_cHOY_103_10483695_M53J6_R1                        5
#  7 CT_cHOY_122_10483695_M53J6_R1                        5
#  8 CT_cHOY_124_10483695_M53J6_R1                        7
#  9 CT_cHOY_125_10483695_M53J6_R1                        5
# 10 CT_cHOY_3_10483695_M53J6_R1                          2
# 11 CT_cHOY_41_10483695_M53J6_R1                         2
# 12 CT_cHOY_42_10483695_M53J6_R1                         2


# Some "good" reads but a lot of single primer for some sets
print(
  OK_sub %>% filter(primer_set=="1325596_93"),width = Inf) 

# # A tibble: 5 × 13
#   primer_set total_primers_sum total_reads only_forward only_reverse
#   <chr>                  <int>       <int>        <int>        <int>
# 1 1325596_93            498191       76245         2833            0
# 2 1325596_93            498191       81085         3055            0
# 3 1325596_93            498191       65229         2432            0
# 4 1325596_93            498191       66662         2518            0
# 5 1325596_93            498191       58520         2406            0
#   only_forward_rc only_reverse_rc forward_reverse forward_rc_reverse
#             <int>           <int>           <int>              <int>
# 1               0               6               0                  0
# 2               0               7               0                  0
# 3               0               7               0                  0
# 4               0               4               0                  0
# 5               0               7               0                  0
#   forward_reverse_rc forward_rc_reverse_rc total_primers
#                <int>                 <int>         <int>
# 1                  3                     0          2842
# 2                  2                     0          3064
# 3                  3                     0          2442
# 4                  2                     0          2524
# 5                  2                     0          2415
#   individuals
#   <chr>
# 1 CT_Tow1_noSQP_455_2_10483695_M53J6_R1
# 2 CT_Tow1_noSQP_455_9_10483695_M53J6_R1
# 3 CT_cART_466_46_10483695_M53J6_R1
# 4 CT_cHOY_8_10483695_M53J6_R1
# 5 CT_Tow1_noSQP_300_12_10483695_M53J6_R1

# primer set with most count
filtered_table %>% filter(total_primers = max(total_primers))
print(filtered_table %>% filter(primer_set == "101680_103"),width = Inf)


# Following Amanda
# 1- Plot the sum over all ind (not mean) for reads with only 1 primer (AmpliconReadCounter seems to add)
#   need to order with higher first - maybe add labelling of this locus on fig (only ones selected)
# 2- Plot with a single individual to see if same loci looking bad
#   need to order with higher first - maybe add labelling of this locus on fig (only ones selected)

sum_single <- count_table %>%
  group_by(primer_set) %>%
  summarise(total_primers_single = sum( c_across(contains("only_")) )) %>%
  arrange(desc(total_primers_single)) %>%
  mutate(primer_set=factor(primer_set, levels=primer_set)) #update the factor levels

limit_read = 1000000
nrow( sum_single %>% filter(total_primers_single>limit_read))

single_count_plot <- ggplot(sum_single, aes(x = primer_set, y = total_primers_single)) +
  geom_point() +
  labs(title = "Primer Set - Single Primer Reads Count (Individuals Sum)",
       x = "Primer Sets",
       y = "Total Single Primers") +
  geom_hline(yintercept = limit_read, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
ggsave(single_count_plot, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/single_count_plot.png")


# Calculate the cumulative sum
df_cumsum <- sum_single %>%
  arrange(total_primers_single) %>% # Arrange by count in ascending order
  mutate(primer_set=factor(primer_set, levels=primer_set)) %>%#update the factor levels
  mutate(
    total_single = sum(total_primers_single), # Calculate the total sales (same for all rows)
    cumulative_count = cumsum(total_primers_single), # Calculate the running total
    cumulative_proportion = cumulative_count / total_single # Calculate the proportion
    # Alternatively, you can combine the last two steps:
    # cumulative_proportion = cumsum(sales_count) / sum(sales_count)
  )
#cum_limit = 0.2
cum_limit = 0.1
nrow( df_cumsum %>% filter(cumulative_proportion>cum_limit))

cumsum_plot <- ggplot(df_cumsum, aes(x = primer_set, y = cumulative_proportion)) +
  geom_point() +
  labs(title = "Primer Set - Single Primer Reads Cumulative Proportion (all ind)",
       x = "Primer Sets",
       y = "Cumulative Proportion") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
ggsave(cumsum_plot, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/cumsum_plot_plot.png")



single_ind <- count_table %>%
  filter(individuals == "CT_cHOY_124_10483695_M53J6_R1") %>%
  group_by(primer_set) %>%
  summarise(total_primers_single = sum( c_across(contains("only_")) )) %>%
  arrange(desc(total_primers_single)) %>%
  mutate(primer_set=factor(primer_set, levels=primer_set)) #update the factor levels

limit_read_ind = 5000
nrow( single_ind %>% filter(total_primers_single>limit_read))

single_count_ind_plot <- ggplot(single_ind, aes(x = primer_set, y = total_primers_single)) +
  geom_point() +
  labs(title = "Primer Set - Single Primer Reads Count (CT_cHOY_124_10483695_M53J6_R1)",
       x = "Primer Sets",
       y = "Total Single Primers") +
  geom_hline(yintercept = limit_read_ind, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10)) 

ggsave(single_count_ind_plot, filename="/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/single_count_ind_plot.png")


# Saving the "to_keep_for_run2", "to_remove_for_run2" list
to_keep <- df_cumsum %>% 
  filter(cumulative_proportion<=cum_limit) %>%
  #sum_single %>% 
  #filter(total_primers_single<=limit_read) %>%
  select(primer_set) %>%
  mutate(
    Name_Forward = paste(primer_set, "_f1", sep = ""),
    Name_Reverse = paste(primer_set, "_r1", sep = "")
    ) 
to_keep <- to_keep %>%  mutate(rownames=rownames(to_keep), PCR_reaction=paste("PCR_",rownames,sep="")) %>%
  select(-rownames)


to_remove <- df_cumsum %>% 
  filter(cumulative_proportion>cum_limit) %>%
  #sum_single %>% 
  #filter(total_primers_single>limit_read) %>%
  select(primer_set) %>%
  mutate(
    Name_Forward = paste(primer_set, "_f1", sep = ""),
    Name_Reverse = paste(primer_set, "_r1", sep = "")
    )

  # Making sure I got the right primer names 
  # cisco_100wind_relax_primers_ordering_format_all_dimers.csv
# To do maybe using the best individual match primers with deltaG but again mean delta for 1 set could be a very bad interaction 
primers_run1 <- fread("/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/tm57/cisco_100wind_relax_primers_ordering_format_all_dimers.csv", header = TRUE, sep = ",")

test1 <- primers_run1 %>% filter(Name %in% to_keep$Name_Forward | Name %in% to_keep$Name_Reverse)
test2 <- primers_run1 %>% filter( ( ! Name %in% to_keep$Name_Forward) & (! Name %in% to_keep$Name_Reverse))

# Save 
out_dir <- "/home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/"
output_filename_keep <- paste0(out_dir,"primers_to_keep_for_run2.csv")
output_filename_remove <- paste0(out_dir,"primers_to_remove_for_run2.csv")

write.csv(to_keep, file = output_filename_keep, row.names = FALSE, quote = FALSE)
write.csv(to_remove, file = output_filename_remove, row.names = FALSE, quote = FALSE)
