#!/usr/bin/env Rscript

# TODO: Incorporate this script into one additional QC report script per project.

# Load required libraries
pacman::p_load(argparse, dplyr, tidyverse, CRISPRcleanR, writexl)

# ------------------------------
# Parse command line arguments
# ------------------------------
parser <- ArgumentParser(
  description = "run_cnv_impact.R: perform a comparison between MAGECK RRA CNV corrected and uncorrected results."
)

parser$add_argument(
  "--mageck-gene-summ-corr",
  help = "Path to MAGECK RRA gene summary for CNV corrected data.",
  type = "character"
)

parser$add_argument(
  "--mageck-gene-summ-uncorr",
  help = "Path to MAGECK RRA gene summary for CNV uncorrected data.",
  type = "character"
)

parser$add_argument(
  "--fdr-threshold",
  help = "FDR rate threshold at which genes are called as significantly exerting a loss/gain-of-fitness effect.",
  type = "double"
)

parser$add_argument(
  "--label",
  help = "Experiment name",
  type = "character"
)

parser$add_argument(
  "--out-cnv-impact-stats",
  help = "Path to save CNV impact statistics (Excel file)",
  type = "character"
)

parser$add_argument(
  "--out-cnv-impact-plots",
  help = "Path to save CNV impact plots (PDF)",
  type = "character"
)

args <- parser$parse_args()

# ------------------------------
# Extract arguments
# ------------------------------
uncorr_res <- args$mageck_gene_summ_uncorr
corr_res <- args$mageck_gene_summ_corr
fdr_thres <- args$fdr_threshold
label <- args$label
out_stat <- args$out_cnv_impact_stats
out_plot <- args$out_cnv_impact_plots

# ------------------------------
# Validate input files
# ------------------------------
if (!file.exists(uncorr_res)) {
  stop("MAGECK uncorrected file does not exist.")
}
if (!file.exists(corr_res)) {
  stop("MAGECK corrected file does not exist.")
}

cat("[INFO] Running CNV impact analysis for", label, "\n")


# ------------------------------
# Create PDF for plots
# ------------------------------
pdf(out_plot, width = 6, height = 6)


# ------------------------------
# Assess CNV correction impact
# ------------------------------
cnv_corr_impact <- ccr.impactOnPhenotype(
  MO_uncorrectedFile = uncorr_res,
  MO_correctedFile   = corr_res,
  sigFDR             = fdr_thres,
  expName            = label,
  display            = FALSE
)

dev.off()

cat("[INFO] CNV impact analysis for", label, "finished: STATUS OK\n")


# ------------------------------
# Extract numeric outputs
# ------------------------------
numeric_outputs <- cnv_corr_impact[sapply(cnv_corr_impact, function(x) is.numeric(x) && length(x) == 1)]

# General statistics from CNV corrected and uncorrected outputs
res_stats <- enframe(unlist(numeric_outputs))

# Contingency table from depletion and enrichment profiles after CNV correction
gene_counts <- cnv_corr_impact$geneCounts

# Genes whose fitness effect has been distorted by CNV correction
distorsion <- as.data.frame(cnv_corr_impact$distorsion)

# Genes whose fitness effect has been attenuated by CNV correction
impact <- as.data.frame(cnv_corr_impact$impact)


# ---------------------------------------
# Check each output object before saving
# ---------------------------------------
for (obj_name in c("res_stats", "gene_counts", "distorsion", "impact")) {
  obj <- get(obj_name)

  if (is.null(obj)) {
    cat("[WARN]", obj_name, "no tiene resultados disponibles (NULL)\n")
    obj <- data.frame(Mensaje = paste("No hay resultados para", obj_name))
  } else if (!is.data.frame(obj)) {
    obj <- tryCatch(as.data.frame(obj), error = function(e) {
      cat("[WARN]", obj_name, "no pudo convertirse a data.frame\n")
      data.frame(Mensaje = paste("No hay resultados válidos para", obj_name))
    })
  }

  assign(obj_name, obj)
}

# ------------------------------
# Collect all outputs
# ------------------------------
cnv_stats <- list(
  CNV_correction_impact_stats    = res_stats,
  CNV_correction_gene_counts     = gene_counts,
  CNV_correction_distorsion      = distorsion,
  CNV_correction_impact          = impact
)

# ------------------------------
# Save results to Excel
# ------------------------------
write_xlsx(cnv_stats, path = out_stat)
