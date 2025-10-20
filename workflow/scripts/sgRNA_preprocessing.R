#!/usr/bin/env Rscript

# Libraries
pacman::p_load(argparse, dplyr, tidyverse, glue)

# Parse data
parser <- ArgumentParser(description = "sgRNA_preprocessing perform identification, filtering and reporting of low-abundant sgRNA.")
parser$add_argument("--input", "-i", help = "sgRNA raw counts file (TSV)", type = "character")
parser$add_argument("--out-count-filt", help = "sgRNA filtered counts (TSV)", type = "character")
parser$add_argument("--out-count-filt-NCT", help = "sgRNA filtered counts with off-target sgRNAs (TSV)", type = "character")
parser$add_argument("--out-filter-report", help = "Report ", type = "character")
parser$add_argument("--out-sgrna-removed", help = "Full description of sgRNA removed.", type = "character")
parser$add_argument("--control-samples", help = "Control samples in the experimental design.", type = "character")
parser$add_argument("--treat-samples", help = "Trteated/Second time samples in the experimental design.", type = "character")
parser$add_argument("--project", help = "Label to use in results. Eg: Project name.", type = "character")
parser$add_argument("--extra", help = "Extra parameters", type = "character")
args <- parser$parse_args()

# Extract args
input_file <- args$input
output_count_filt <- args$out_count_filt
output_count_filt_NCT <- args$out_count_filt_NCT
output_filter_report <- args$out_filter_report
output_sgrna_removed <- args$out_sgrna_removed
controls <- args$control_samples
treated <- args$treat_samples
project <- args$project
extra <- args$extra

message(glue("Input file: {input_file}"))

controls <- unlist(strsplit(controls, ","))
treated <- unlist(strsplit(treated, ","))
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
zero_counts <- sgRNA_counts %>% filter(rowSums(select(., -sgRNA, -Gene), na.rm = TRUE) == 0)

# 3.2 Identify sgRNAs with zero counts in control samples and any counts in treated samples.
zero_counts_ctrl <- sgRNA_counts %>%
  filter(rowSums(select(., all_of(controls)), na.rm = T) == 0 & rowSums(select(., all_of(treated)), na.rm = TRUE) > 0)

# 3.3 Indentify sgRNAs with mean abundance lower than 30 in control samples.
mean_lower_thres_ctrl <- sgRNA_counts %>%
  filter(rowSums(select(., -sgRNA, -Gene), na.rm = TRUE) > 0) %>%
  filter(!(rowSums(select(., all_of(controls)), na.rm = TRUE) == 0 & rowSums(select(., all_of(treated)), na.rm = TRUE) > 0)) %>%
  # TODO: Make the threshold a parameter.
  filter(rowMeans(select(., all_of(controls)), na.rm = TRUE) < 30)


# 3.4 Filter raw counts to remove sgRNAs with zero counts in all samples and in control samples.
raw_counts_filt <- sgRNA_counts %>%
  filter(!sgRNA %in% c(zero_counts$sgRNA, zero_counts_ctrl$sgRNA, mean_lower_thres_ctrl$sgRNA))


# 4. Write general report for sgRNA filtering.
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
), output_filter_report)


# 5. Write full report of lost sgRNAs causes at sgRNA and gene level.

# 5.1 Label the cause of removal for each sgRNA.
sgRNA_counts <- sgRNA_counts %>%
  mutate(lost_cause = case_when(
    sgRNA %in% zero_counts$sgRNA ~ "Zero counts in all samples",
    sgRNA %in% zero_counts_ctrl$sgRNA ~ "Empty in T0/Control samples",
    sgRNA %in% mean_lower_thres_ctrl$sgRNA ~ "Low abundant in T0/Control samples",
    TRUE ~ NA
  ))

# 5.2 Create full report of lost sgRNAs
if (n_zero_all == 0 & n_zero_ctrl == 0 & n_mean_lower_thres_ctrl == 0) {
  message("No sgRNAs were filtered - OK!")
} else {
  filtering_report <- sgRNA_counts %>%
    filter(!is.na(lost_cause)) %>%
    group_by(Gene) %>%
    summarise(
      n_lost_sgRNAs = n(),
      lost_sgRNAs_ids = str_c(sgRNA, collapse = ","),
      lost_sgRNA_cause = str_c(lost_cause, collapse = ",")
    ) %>%
    arrange(desc(n_lost_sgRNAs))
}

# 6. Generate processed library file with and without control sgRNAs.

# 6.1 sgRNA counts without control sgRNAs
sgRNA_counts_processed <- sgRNA_counts %>%
  filter(is.na(lost_cause) & !grepl("Non-Targeting_Control", Gene)) %>%
  select(-lost_cause)

# 6.2 sgRNA counts with control sgRNAs
sgRNA_counts_processed_NCT <- sgRNA_counts %>%
  filter(is.na(lost_cause)) %>%
  select(-lost_cause)


# 7. Export results

# 7.1 Export detailed sgRNA filtering report.
if (exists("filtering_report")) {
  write.table(filtering_report, output_sgrna_removed, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
  write_lines("No sgRNAs were filtered - OK!", output_sgrna_removed)
}

# 7.2 Export filtered raw counts with control sgRNAs.
write.table(sgRNA_counts_processed_NCT, output_count_filt_NCT, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 7.3 Export filtered raw counts without control sgRNAs.
write.table(sgRNA_counts_processed, output_count_filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
