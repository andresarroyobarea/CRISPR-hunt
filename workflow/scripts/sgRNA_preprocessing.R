#!/usr/bin/env Rscript

# Libraries
pacman::p_load(argparse, dplyr, tidyverse, glue)

# Parse data
parser <- ArgumentParser(description= 'sgRNA_preprocessing perform identification, filtering and reporting of low-abundant sgRNA.')
parser$add_argument('--input', '-i', help = 'sgRNA raw counts file (TSV)', type = "character")
parser$add_argument('--out-count-filt', help = 'sgRNA filtered counts (TSV)', type = "character")
parser$add_argument('--out-count-filt-NCT', help = 'sgRNA filtered counts with off-target sgRNAs (TSV)', type = "character")
parser$add_argument('--out-filter-report', help = 'Report ', type = 'character')
parser$add_argument('--out-sgrna-removed', help = 'Full description of sgRNA removed.', type = 'character')
parser$add_argument('--control-samples', help = 'Control samples in the experimental design.', type = "character")
parser$add_argument('--treat-samples', help = 'Trteated/Second time samples in the experimental design.', type = "character")
parser$add_argument('--project', help = 'Label to use in results. Eg: Project name.', type = "character")
parser$add_argument('--extra', help = 'Extra parameters', type = "character")
args <- parser$parse_args()

# Extract args
input_file <- args$input
output_count <- args$out_count_filt
output_count_NCT <- args$out_count_filt_NCT
output_filter_report <- args$out_filter_report
output_sgrna_removed <- args$out_sgrna_removed
controls <- args$control_samples
treated <- args$treat_samples
project <- args$project
extra <- args$extra

message(glue("Input file: {input_file}"))

controls <- unlist(strsplit(c(controls), ","))
treated <- unlist(strsplit(c(treated), ","))
print(controls)
print(treated)

# 1. Read sgRNA raw counts file
if (!file.exists(input_file)) {
  stop(paste0("Input file does not exist: ", input_file))
} else {
  message("Reading input file: ", input_file)
}

sgRNA_counts <- read.delim(input_file, sep = "\t", header = TRUE)

# 2. Check if the first two columns are 'sgRNA' and 'Gene'.
if (!identical(colnames(sgRNA_counts)[1:2], c("sgRNA", "Gene"))) {
  stop(paste0("Count file must have 'sgRNA' and 'Gene' as the first two columns."))
}

# 3. Remove sgRNAs with zero counts in all samples.
# 3.1 Identify sgRNAs with zero counts in all samples.
zero_counts <- sgRNA_counts %>% filter(rowSums(select(., -sgRNA, -Gene)) == 0)

# 3.2 Identify sgRNAs with zero counts in control samples and any counts in treated samples.
zero_counts_ctrl <- sgRNA_counts %>% 
    filter(rowSums(select(., all_of(controls))) == 0 & rowSums(select(., all_of(treated))) > 0)

# 3.3 Indentify sgRNAs with mean abundance lower than 30 in control samples.
mean_lower_thres_ctrl <- sgRNA_counts %>%
    filter(rowSums(select(., -sgRNA, -Gene)) > 0) %>%
    filter(!(rowSums(select(., all_of(controls))) == 0 & rowSums(select(., all_of(treated))) > 0)) %>%
    # TODO: Make the threshold a parameter.
    filter(rowMeans(select(., all_of(controls))) < 30)


# 3.4 Filter raw counts to remove sgRNAs with zero counts in all samples and in control samples.
raw_counts_filt <- sgRNA_counts %>% 
    filter(!sgRNA %in% c(zero_counts$sgRNA, zero_counts_ctrl$sgRNA, mean_lower_thres_ctrl$sgRNA))

# 4. Write report of sgRNA filtering.
total_sgrnas <- nrow(sgRNA_counts)
n_zero_all <- nrow(zero_counts)
n_zero_ctrl <- nrow(zero_counts_ctrl)
n_mean_lower_thres_ctrl <- nrow(mean_lower_thres_ctrl)
n_after_filter <- nrow(raw_counts_filt)

write_lines(glue(
  "sgRNA Filtering Report\n",
  "----------------------\n",
  "Total sgRNAs in library: {total_sgrnas}\n",
  "sgRNAs with 0 counts in all samples: {n_zero_all} ({round(n_zero_all / total_sgrnas * 100, 2)}%)\n",
  "sgRNAs with zero counts in all T0/Control samples but non-zero counts in any later sample: {n_zero_ctrl} ({round(n_zero_ctrl / total_sgrnas * 100, 2)}%)\n",
  "sgRNAs with low abundance in T0/Control samples (mean abundance < 30): {n_mean_lower_thres_ctrl} ({round(n_mean_lower_thres_ctrl / total_sgrnas * 100, 2)}%)\n",
  "sgRNAs after filtering: {n_after_filter} ({round(n_after_filter / total_sgrnas * 100, 2)}%)\n"
), filt_sum)

# NEW ---
sgRNA_counts %>%
  mutate(cause = case_when(
    sgRNA %in% zero_counts$sgRNA ~ "Zero counts in all samples",
    sgRNA %in% zero_counts_ctrl$sgRNA ~ "Empty in T0/Control samples",
    sgRNA %in% mean_lower_thres_ctrl$sgRNA ~ "Low abundant in T0/Control samples",
    TRUE ~ NA
  ))
# ---


# Processed version

# NCT version.




# 5. Write report of lost sgRNAs at the gene level.
if (nrow(zero_counts) == 0) {
  message("No sgRNAs with zero counts in all samples - OK!")
} else {
  zero_counts_all_by_gene <- zero_counts %>% 
    group_by(Gene) %>% 
    summarise(lost_sgRNAs = n()) %>% 
    arrange(desc(lost_sgRNAs))
}

# Export gene summary report
write.csv(zero_counts_all_by_gene, gene_summ, sep = "\t", row.names = FALSE, quote = FALSE)

# Export filtered raw counts.
write.csv(raw_counts_filt, output_count, sep = "\t", rownames = FALSE, quote = FALSE)
