#!/usr/bin/env Rscript

# Libraries
pacman::p_load(argparse, dplyr, tidyverse, tibble, CRISPRcleanR)

# Utility functions
plot_qc_curves <- function(feature_lfc, label) {

    """
    Function to plot ROC and Precision/Recall 
    curves from CRISPRcleanR results
    """

    # ROC curve
    ccr.ROC_Curve(
        FCsprofile = feature_lfc,
        positives = bagel2_essentials,
        negatives = bagel2_non_essentials,
        display = FALSE,
        FDRth = 0.05,
        expName = label)

    # Precision/Recall curve
    ccr.PrRc_Curve(
        FCsprofile = feature_lfc,
        positives = bagel2_essentials,
        negatives = bagel2_non_essentials,
        display = FALSE,
        FDRth = 0.05, # Turn this parameter into an option?
        expName = label)

}

# Parse data
# Parse data
parser <- ArgumentParser(description= 'run_crisprcleanR_QC perform QC in CRISPRcleanR CNV corrected data.')
parser$add_argument('--sgrna-summ', help = 'MAGECK sgRNA summary (RRA/MLE) (TSV)', type = "character")
parser$add_argument('--gene-summ', help = 'MAGECK sgRNA summary (RRA/MLE) (TSV)', type = "character")
parser$add_argument('--lib-type', help = 'Custom or CRIPRcleanR default library.', type = 'character')
parser$add_argument('--sgrna-library', help = 'sgRNA library file: "custom" or "built-in" (CRISPRcleanR)', type = 'character')
parser$add_argument('--sel-type', help = 'CRISPRscreening selection type: positive or negative.', type = 'character')
parser$add_argument('--essen-genes', help = 'Essential genes (default: BAGEL2 essential genes)', type = 'character')
parser$add_argument('--sgrna-library', help = 'Non-essential genes (default: BAGEL2 non-essential genes)', type = 'character')
parser$add_argument('--label', help = 'Project name.', type = 'character')

args <- parser$parse_args()

# Extract args
sgrna_summ <- args$sgrna_summary
gene_summ <- args$gene_summary
output_file <- args$qc_report
lib_type <- args$library_type
sgrna_lib <- args$sgrna_library
selection_type <- args$selection_type
essential_genes <- args$essential_genes
non_essential_genes <- args$non_essential_genes
label <- args$label

# CRISPRcleanR custom data
# Essential gene sets
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

# Chequear si hay versiones nuevas y más actualizadas de estos genes
core_fitness_genes <- list(
    Ribosomal_Proteins=EssGenes.ribosomalProteins,
    DNA_Replication = EssGenes.DNA_REPLICATION_cons,
    RNA_polymerase = EssGenes.KEGG_rna_polymerase,
    Proteasome = EssGenes.PROTEASOME_cons,
    Spliceosome = EssGenes.SPLICEOSOME_cons,
    CFE = bagel2_essentials,
    non_essential = bagel2_non_essentials)

# Import sgRNA and gene summaries.
sgrna_input <- read.delim(sgrna_summ, sep = "\t", header = TRUE, row.names = FALSE)
gene_input <- read.delim(gene_summ, sep = "\t", header = TRUE, row.names = FALSE)
essential_genes <- read.delim(essential_genes, sep = "\t", header = TRUE, row.names = FALSE)
non_essential_genes <- read.delim(non_essential_genes, sep = "\t", header = TRUE, row.names = FALSE)

# Import library file or load CRISPRcleanR library.
# Format Checkings are not repeated here because they were done in run_crisprcleanR.R
if (lib_type == "custom") {

    # Load user-provider custom library
    message("Using custom library": sgrna_lib)

    # Import KY_Library_v1.0 for further checkings.
    ky_library <- get(data(KY_Library_v1.0))

    # TODO: Chequear si tenemos los sgRNAs ids como rownames y como columna CODE.
    lib_file = read.csv(sgrna_lib, sep = "\t", header = TRUE) 
    
} else {

    # Load built-in library
    message("Using internal CRISPRcleanR library: ", sgrna_lib)

    lib_file <- get(data(sgrna_lib))

}

# ---- Preprocessing ---- #
sgrna_lfc <- sgrna_input %>% select(sgRNA, log2FC) %>% deframe()


if (selection_type == "negative") {

    # Negative lfc is selected
    gene_lfc <- gene_input %>% select(id, neg.lfc) %>% deframe()

} else {
    
    # Positive lfc is selected
    gene_lfc <- gene_input %>% select(id, pos.lfc) %>% deframe()

}


# BAGEL2 essential genes.
bagel2_essentials <- essential_genes %>% select(GENE) %>% pull()
bagel2_essentials <- ccr.genes2sgRNAs(lib_file, bagel2_essentials)

# BAGEL2 non-essential genes.
bagel2_non_essentials <- non_essential_genes %>% select(GENE) $>% pull()
bagel2_non_essentials <- ccr.genes2sgRNAs(lib_file, bagel2_non_essentials)


# ---- Run CRISPRcleanR QC ---- #
pdf(output_file)

if (exits("sgrna_lfc")) {

    # Run ROC and Precision/Recall curves if sgRNA LFC file exists
    plot_qc_curves(sgrna_lfc, paste0(label, " (sgRNA-level)"))
}

if (exists("gene_lfc")) {

    # Run ROC and Precision/Recall curves if sgRNA LFC file exists
    plot_qc_curves(gene_lfc, paste0(label, " (gene-level)"))

    # Depletion profile at gene level
    ccr.VisDepAndSig(
        FCsprofile = gene_lfc,
        SIGNATURES = core_fitness_genes,
        TITLE = label,
        pIs = 6,
        nIs = 7,
        plotFCprofile = TRUE)
}

dev.off()




# AQUI SE DEBE AÑADIR UN BUCLE QUE SE EJECUTE SI HAY INFROMACIÓN GENOMICA DE LA LINEA CELULAR EN 
# CUESTION CON ccr.perf_statTests, ccr.perf_distributions y ccr.RecallCurves


# LUEGO DEBE HABER OTRO BUCLE QUE SE EJECUTE PARA COMPARAR RESULTADOS DE MAGECK CORREGIDOS Y SIN CORREGIR.
# OBJETIVO: VER EL EFECTO DE LA CORRECION.


