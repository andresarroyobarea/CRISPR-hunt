#!/usr/bin/env Rscript

# =========================
# Libraries
# =========================
pacman::p_load(argparse, dplyr, tidyverse, purrr, CRISPRcleanR)

# =========================
# Argument parser
# =========================
parser <- ArgumentParser(description = "run_crisprcleanR performs unsupervised CNV correction in CRISPR screening data.")

parser$add_argument("--input", "-i", help = "sgRNA raw counts file (TSV)", type = "character")
parser$add_argument("--out-count", "-o", help = "Output sgRNA CNV-corrected and normalized counts file.", type = "character")
parser$add_argument("--out-lfc", help = "Output sgRNA CNV-corrected and normalized log2 fold-change file.", type = "character")
parser$add_argument("--lib-type", help = "Custom or CRISPRcleanR default library.", type = "character")
parser$add_argument("--sgrna-library", help = 'sgRNA library file: "custom" or "built-in" (CRISPRcleanR)', type = "character")
parser$add_argument("--norm-method", help = "Normalization method prior to CNV correction.", type = "character")
parser$add_argument("--exp-design", help = "Experimental design to calculate log2FC for CNV-corrected sgRNA count data.", type = "character")
parser$add_argument("--min-reads", help = "Minimum number of reads per sgRNA to be included.", type = "double")
parser$add_argument("--min-genes", help = "Minimal number of different genes targeted by sgRNAs in a biased segment to perform count correction.", type = "double")
parser$add_argument("--control-samples", help = "Control samples in the experimental design.", type = "character")
parser$add_argument("--treat-samples", help = "Treated/Condition samples in the experimental design.", type = "character")
parser$add_argument("--label", help = "Label to use in results (e.g. project name).", type = "character")
parser$add_argument("--outdir", help = "Folder to save the results.", type = "character")
parser$add_argument("--extra", help = "Extra parameters, e.g., DNAcopy arguments.", type = "character")

args <- parser$parse_args()

# =========================
# Extract arguments
# =========================
input_file <- args$input
output_count <- args$out_count
output_lfc <- args$out_lfc
lib_type <- args$lib_type
sgrna_lib <- args$sgrna_library
norm_method <- args$norm_method
exp_design <- args$exp_design
min_reads <- args$min_reads
min_genes <- args$min_genes
controls <- args$control_samples
treated <- args$treat_samples
label <- args$label
outdir <- args$outdir
extra <- args$extra

# =========================
# Control and treated samples info
# =========================
print(paste("Control samples:", controls))
n_controls <- length(unlist(strsplit(controls, ",")))
print(paste("Number of control samples:", n_controls))

controls <- strsplit(controls, ",")[[1]] %>% trimws()
treated <- strsplit(treated, ",")[[1]] %>% trimws()

# =================================
# Raw counts importing and checking
# =================================
# TODO: CHECK
raw_counts <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

expected_cols <- c("sgRNA", "gene")
actual_cols <- colnames(raw_counts)[1:2]

if (!all(actual_cols == expected_cols)) {
  warning(sprintf(
    "Warning: Column names differ from expected (%s). Automatically renaming columns to 'sgRNA' and 'gene'.",
    paste(actual_cols, collapse = ", ")
  ))
  colnames(raw_counts)[1:2] <- expected_cols
}

# =====================================
# sgRNA library importing and checking
# =====================================
if (lib_type == "built-in") {
  message("Using built-in CRISPRcleanR sgRNA library.")
  sgRNA_library <- data(sgrna_lib)
} else if (lib_type == "custom") {
  # Load user-provided custom library
  message(paste0("Using custom library: ", sgrna_lib))
  sgRNA_library <- read.table(sgrna_lib, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Validate format matches KY_Library_v1.0 structure
  required_cols <- c("CODE", "GENES", "EXONE", "CHRM", "STRAND", "STARTpos", "ENDpos", "seq")

  if (!all(required_cols %in% colnames(sgrna_lib))) {
    stop(paste0(
      "Error: Custom library must follow the KY_Library_v1.0 format and have columns: ",
      paste(required_cols, collapse = ", ")
    ))
  }
} else {
  stop("Error: lib_type must be either 'custom' or 'built-in'.")
}

# =====================================
# Apply CRISPRcleanR CNV correction
# =====================================

# 1. Normalization and logFC calculation.
norm_and_fcs <- ccr.NormfoldChanges(
  Dframe = raw_counts,
  saveToFig = TRUE,
  min_reads = min_reads,
  EXPname = label,
  libraryAnnotation = sgRNA_library,
  n_controls = n_controls,
  display = TRUE,
  method = norm_method
)

print("NORM FOLD CHANGES - OK!!!")
print("First rows of normalized counts:")
print(head(norm_and_fcs[["norm_counts"]]))
print("First rows of log2 fold changes:")
print(head(norm_and_fcs[["logFCs"]]))

# 2. Genomic sorting of sgRNAs’ log fold changes.
gw_sorted_fcs <- ccr.logFCs2chromPos(
  foldchanges = norm_and_fcs[["logFCs"]],
  libraryAnnotation = sgRNA_library
)

print("sgRNAs L2FC GENOMIC SORTING - OK!!!")
print("First 5 rows of mapped log fold changes:")
print(head(gw_sorted_fcs))

# 3. Unsupervised identification and correction of gene independent cell responses to CRISPR-Cas9 targeting
corrected_fcs <- ccr.GWclean(
  gwSortedFCs = gw_sorted_fcs,
  label = label,
  display = TRUE,
  saveTO = TRUE,
  ignoredGenes = NULL,
  min.ngenes = min_genes
)
print("sgRNAs L2FC CORRECTION - OK!!!")

# 4. Correction of sgRNA treatment counts for gene independent responses to CRISPR-Cas9 targeting
corrected_counts <- ccr.correctCounts(
  CL = label,
  normalised_counts = norm_and_fcs[["norm_counts"]],
  correctedFCs_and_segments = corrected_fcs,
  libraryAnnotation = sgRNA_library,
  minTargetedGenes = min_genes,
  OutDir = outdir,
  ncontrols = n_controls,
)

# =====================================================
# CNV-corrected log2 fold changes per treatment sample
# =====================================================
# Prepare CNV corrected log2 fold changes to be used by BAGEL2.

# ---1. Pseudo-count addition to avoid log2(0) ---
pseudocount <- 0.1
corrected_counts <- corrected_counts %>%
  mutate(across(all_of(c(controls, treated)), ~ .x + pseudocount))


# --- Calculate log2 fold changes based on experimental design ---

if (exp_design == "global") {
  log2fc_cnv_corrected <- corrected_counts %>%
    rowwise() %>%
    mutate(ctrl_mean = mean(c_across(all_of(controls)))) %>%
    ungroup() %>%
    mutate(across(all_of(treated), ~ log2(.x / ctrl_mean))) %>%
    select(sgRNA, gene, all_of(treated)) %>%
    as.data.frame()
} else if (exp_design == "paired") {
  if (length(controls) != length(treated) & exp_design == "paired") {
    stop("Error: In paired designs, the number of control and exposed samples must be the same.")
  }

  log2fc_cnv_corrected <- corrected_counts %>%
    select(sgRNA, gene) %>%
    bind_cols(
      map2_dfc(treated, controls, ~ log2(corrected_counts[[.x]] / corrected_counts[[.y]])) %>%
        set_names(treated) # asigna nombres de columnas directamente
    ) %>%
    as.data.frame()
} else {
  stop("Error: exp_design must be either 'global' or 'paired'.")
}

# =====================================================
# Export results
# =====================================================
# CNV-corrected counts
write.table(corrected_counts, file = output_count, sep = "\t", quote = FALSE, row.names = FALSE)

# CNV-corrected log2 fold changes
write.table(log2fc_cnv_corrected, file = output_lfc, sep = "\t", quote = FALSE, row.names = FALSE)

print("CRISPRcleanR CNV correction completed successfully.")
