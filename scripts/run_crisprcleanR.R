#!/usr/bin/env Rscript

# Libraries
pacman::p_load(argparse, dplyr, tidyverse, CRISPRcleanR)

# Parse data


# Checks












# Run CRISPRcleanR

# IMPORTANT: 
# ncontrols: Number. it takes the number of columns to take as controls after gene.  
# A. Controls must be always after the sgRNA and gene columns for CRISPRcleanR.
# B. CRISPRcleanR takes the mean of normalized control sgRNA counts to calculate log2FC.
# Eg: control1, control2, treat1, treat2
# norm_sgrna_1: 20, 30, 18, 40
# lFC sgRNA1 treat1 = log2(18 / mean(c(20, 30)))
# lFC sgRNA1 treat2 = log2(40 / mean(c(20, 30)))
# Conclusion: First, provide controls after sgRNA and gene columns. Second, specify the number of
# control samples.

# 1. Normalisation of sgRNA and fold change computation.
norm_and_logfcs <- ccr.NormfoldChanges(
    filename = , 
    libraryAnnotation = , 
    min_reads =  , 
    EXPname = , 
    saveToFig = TRUE, 
    outdir =  ,
    ncontrols = , 
    method = )


# 2. Genomic sorting of sgRNA's log fold changes
mapped_logfcs <- ccr.logFCs2chromPos(
    foldchanges = norm_and_logfcs[["logFCs"]],
    libraryAnnotation = 
)

# Este output se podría guardar también.

# 3. Unsupervised identification and correction of gene independent cell responses to CRISPR-Cas9 targeting.
corrected_fcs <- ccr.GWclean(
    gwSortedFCs = mapped_logfcs,
    label = # Project name: Parameter,
    saveTO = # Outdir for results,
    ignoredGenes = # BAGEL2 known essential genes? RIbosomal?
    min.ngenes = 	
)

# El resto de parametros de la función previa se dejan por defecto.
# Parece que los log2FC corregidos y los segmentos no se guardan por defecto ---> Pensar si sacarlos como outputs.

# 4. Correction of sgRNA counts for gene independent responses to CRISPR-Cas9 targeting.
corrected_counts <- ccr.correctCounts(
    CL = , # Nombre del proyecto,
    normalised_counts = norm_and_logfcs[["norm_counts"]],
    correctedFCs_and_segments = corrected_fcs,
    libraryAnnotation = ,
    minTargetedGenes = # Parameter min_genes,
    OutDir = # parameter outdir,
    ncontrols = # Parameter n_controls
)

write.csv(corrected_counts, , sep = "\t", rownames = FALSE, quote = FALSE)




# Guardar fichero final para mapear la regla de Snakemake



