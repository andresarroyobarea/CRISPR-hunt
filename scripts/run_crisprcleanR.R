#!/usr/bin/env Rscript

# Libraries
pacman::p_load(argparse, dplyr, tidyverse, CRISPRcleanR)

# Parse data
parser <- ArgumentParser(description= 'run_crisprcleanR perform unsupervised CNV correction in CRISPR screening data.')
parser$add_argument('--input', '-i', help = 'sgRNA raw counts file (TSV)', type = "character")
parser$add_argument('--output', '-o', help = 'Output sgRNA CNVs corrected and normalized counts file.', type = "character")
parser$add_argument('--lib-type', help = 'Custom or CRIPRcleanR default library.', type = 'character')
parser$add_argument('--sgrna-library', help = 'sgRNA library file: "custom" or "built-in" (CRISPRcleanR)', type = 'character')
parser$add_argument('--norm-method', help = 'Normalizacion method prior CNV correction.', type = "character")
parser$add_argument('--min-reads', help = 'Minimal number of sgRNA counts accross all samples to retain the feature', type = "double")
parser$add_argument('--min-genes', help = 'Minimal number of different genes targeted by sgRNAs in a biased segment to perform count correction.', type = "double")
parser$add_argument('--n-control-samples', help = 'Number of controls replicates to be considered in CRISPRcleanR calculations.', type = "double")
parser$add_argument('--label', help = 'Label to use in results. Eg: Project name.', type = "character")
parser$add_argument('--outdir', help = 'Folder to save the results.', type = "character")
parser$add_argument('--extra', help = 'Extra parameters. Eg: DNAcopy arguments', type = "character")

args <- parser$parse_args()

# Extract args
input_file <- args$input
norm_cnv_counts <- args$output
lib_type <- args$library_type
sgrna_lib <- args$sgrna_library
norm_method <- args$norm_method
min_reads <- args$min_reads
min_genes <- args$min_genes
n_control <- args$n_controls
label <- args$project
outdir <- args$outdir
extra <- args$extra


# Import raw counts.
raw_counts <- read.delim(input_file, sep = "\t", header = TRUE, row.names = FALSE)

if (c("sgRNA", "gene") != colnames(raw_counts)[1:2]){
    stop(paste0("Count file must have 'sgRNA' and 'gene' as the first two columns."))
}

# Import library file or load CRISPRcleanR library
# 1. CHECK LIBRARY FILE FORMAT.
if (lib_type == "custom") {

    # Load user-provider custom library
    message("Using custom library": sgrna_lib)

    # Import KY_Library_v1.0 for further checkings.
    ky_library <- get(data(KY_Library_v1.0))

    # TODO: Chequear si tenemos los sgRNAs ids como rownames y como columna CODE.
    lib_file = read.csv(sgrna_lib, sep = "\t", header = TRUE) 

    # Validate format matches KY_Library_v1.0.
    required_cols = c("CODE", "GENES", "EXONE", "CHRM", "STRAND", "STARTpos", "ENDpos", "seq")
    if (all(required_cols != colnames(ky_library)){
        stop(paste0(
            "Error: Custom library must follow the KY_Library_v1.0 format and have columns:" required_cols,
        ))
    }
    
} else {

    # Load built-in library
    message("Using internal CRISPRcleanR library: ", sgrna_lib)

    lib_file <- get(data(sgrna_lib))

}

# 2. CHECK sgRNA duplicates.
dup_sgrna_ids <- lib_file$CODE[duplicated(lib_file$CODE)]

if length(dup_sgrna_ids) > 0 {

    # Filter sgRNA IDs duplicated records
    dup_rows <- s2 %>% filter(CODE %in% dup_ids)

    # Verify if sgRNA sequence is similar in duplicated cases (real duplicates) or
    # only the ID is duplicate (false duplicate).
    # sgRNA with same ID.
    incons_sgrna_ids <- dup_rows %>%
        group_by(CODE) %>%
        summarise(n_unique_seq = n_distinct(seq)) %>%
        filter(n_unique_seq > 1) %>%
        pull(CODE)
    
    # Report sgRNA with the same ID and different sequence to check them
    if (length(incons_sgrna_ids) > 0) {
        stop(paste0(
            "Error: Next sgRNA IDs are duplicated with different sgRNA sequence:",
            paste(incons_sgrna_ids, collapse = ","), ".", "Please, check your library sgRNA IDs.")
            )
    } else {
        message("Warning: There are duplicated sgRNA IDs with the same sequence. The first record was preserved.")

        lib_file <- lib_file[!duplicated(lib_file$CODE), ]
    }

    

} else {
    message("There were not duplicate sgRNA IDs - OK!")
}



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



