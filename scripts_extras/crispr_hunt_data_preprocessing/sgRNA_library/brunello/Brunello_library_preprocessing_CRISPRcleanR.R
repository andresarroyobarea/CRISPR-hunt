# --------------------------------------------------------------------------------------------------------#
# Script: Brunello_library_preprocessing.R
# Description: Preprocessing raw Brunello Library Target Genes to match CRISPRcleanR required input format.
# Author: Andrés Arroyo Barea
# Date: 2025-07-01
# --------------------------------------------------------------------------------------------------------#
pacman::p_load(tidyverse, tibble, here, janitor, HGNChelper, biomaRt, msigdbr, purrr, dplyr)

# -------------------- #
# ----- Objective ---- #
# -------------------- #

# To use a custom library with CRISPRcleanR, it must follow the format used in the KY_Library_v1.0 library
# This library and the 8 variables that should be included are here: https://github.com/francescojm/CRISPRcleanR/blob/master/Reference_Manual.pdf (page 85).

# These variables are:
# CODE alphanumerical identifier of the sgRNAs (string).
# GENES: vector containing the HGNC symbols of the genes targeted by the sgRNA under consideration (string).
# EXONE: vector containing the gene exon targeted by the sgRNA under consideration (these should include the prefix "ex" followed by the number (string).
# CHRM: vector the chromosome of the gene targeted by the sgRNA under consideration (X and Y chromosome should be specified as "X" and "Y" (string).
# STRAND: vector containing the strand targeted by the sgRNA under consideration ("+" or "-") (string).
# STARTpos: vector containing the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration (numeric).
# ENDpos: vector containing the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration (numeric).
# seq: sgRNA sequence (string)

# -------------------- #
# --- 1. Raw Data ---- ####
# -------------------- #

# Brunello Library Target Genes data was downloaded from https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt
# Data was downloaded 2025-07-01.
# A second download was done on 2025-09-24 to read properly Rule 2 Set Score values because they were not well formatted in the previous version.
brunello_raw <- read.delim(here("sgRNA_library", "brunello", "raw", "brunello_library_raw_v2.txt"), sep = "\t")
colnames(brunello_raw) <- make_clean_names(colnames(brunello_raw))

# ----------------------------------------------- #
# ---- 2. Raw Brunello Library preprocessing ---- ####
# ----------------------------------------------- #

# OBJECTIVE: The Brunello Library file was uploaded in the ADDgene webpage some years ago and some genes symbols has been 
# updated between this uploaded date and the current date (2025-09-23). They represent ~ 250-300 target gene symbol in the raw
# library (1% total genes in the library.). In spite of it is not a great amount, they were updated to current gene symbols to avoid
# lose potential hit and to have updated information in pathway analysis.

# Check the number of outdate genes. The updated date at the moment of this checking wass the next: "Maps last updated on: Sat Nov 16 10:35:32 2024"
gene_symbol_check <- checkGeneSymbols(
  unique(brunello_raw$target_gene_symbol),
  unmapped.as.na = TRUE,
  species = "human",
  expand.ambiguous = FALSE)

# 17887 genes with updated gene symbol.
gene_symbol_check %>% filter(Approved == TRUE) %>% dplyr::select(x) %>% pull() %>% unique() %>% length()

# 1228 genes without approved gene symbol
gene_symbol_check %>% filter(Approved == FALSE) %>% dplyr::select(x) %>% pull() %>% unique() %>% length()

# 92/1228 has not a suggested new symbol. Mainly, they are lncRNAs without gene symbol. Control sgRNA excluded.
gene_symbol_check %>% filter(Approved == FALSE & is.na(Suggested.Symbol)) %>% dplyr::select(x) %>% pull() %>% unique() %>% length()

# 1135 / 1228 has outdated gene symbol with suggested new gene symbol
gene_symbol_check %>% filter(Approved == FALSE & !is.na(Suggested.Symbol)) %>% dplyr::select(x) %>% pull() %>% unique() %>% length()

# Rules to update gene symbols.
# 1. If gene symbol is approved --> Keep original Brunello Gene Symbol.
# 2. If gene symbol is not approved and suggested symbol is not found --> Keep original Brunello Gene Symbol
# 3. If gene symbol is not approved and suggested symbol is provided --> Update gene symbol.

# ADDITIONAL 1: Explore genes whose symbol is not approved but they dont have a suggested symbol.
not_approved_symbol <- gene_symbol_check %>% filter(Approved == FALSE & is.na(Suggested.Symbol) & !grepl("Control", x)) %>% dplyr::select(x) %>% pull() %>% unique()

# First, genes with provisional identifiers were checked to explore if they have adquire a gene symbol between the brunello library
# last updated and the current state of gene symbols HGNC.

# Conectar a Ensembl
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

# Select LOC genes
loc_genes <- gsub("-.*", "", not_approved_symbol[grepl("LOC", not_approved_symbol)])

# Extract current information about these loc genes.
res_loc <- getBM(
  attributes=c("entrezgene_id", "ensembl_gene_id","hgnc_symbol","description"),
  filters="entrezgene_id",
  values=as.numeric(sub("LOC","",loc_genes)),
  mart=ensembl
)

# LOC genes which no information when last Brunello Library was updated in AddGene and its
# gene symbol ID updated after this moment
loc_to_update <- res_loc %>% 
  filter(hgnc_symbol != "") %>%
  dplyr::select(entrezgene_id, hgnc_symbol) %>%
  mutate(entrezgene_id = paste0("LOC", entrezgene_id)) %>%
  deframe()


# LOC genes without gene symbol in 2025-09-23
res_loc %>% 
  filter(hgnc_symbol == "") %>%
  dplyr::select(entrezgene_id, hgnc_symbol) %>%
  mutate(entrezgene_id = paste0("LOC", entrezgene_id)) %>%
  deframe()
  
# EXtract no-LOC genes. They were added below in the map_genes_vectot.
not_approved_symbol[!grepl("LOC", not_approved_symbol)]

# ADDITIONAL 2: In a first checking, 14 genes symbol were mapped to more than one updated gene symbol. In these cases,
# the sgRNA sequencewere mapped with the human genome+transcriptome to find the best hits (lower E-value and higher 
# identity) and Gene id number was used to identify the right symbol
# OLD --> UPDATED
map_genes_to_update <- c(
  "GIF" = "CBLIF",
  "LOR" = "LORICRIN",
  "SEPT2" = "SEPTIN2",
  "NOV" = "CCN3",
  "QARS" = "QARS1",
  "SARS" = "SARS1",
  "TAZ" = "TAFAZZIN",
  "DEC1" = "DELEC1",
  "MPP6" = "PALS2",
  "CSRP2BP" = "KAT14",
  "AGPAT9" = "GPAT3",
  "MUM1" = "PWWP3A",
  "DUSP27" = "STYXL2",
  "STRA13" = "CENPX",
  # This one is directly in the name
  "LOC400927" = "TPTEP2-CSNK1E",
  # Not found by HGNChelper and no-LOC genes but found in NCBI (MAYBE CHANGES AFTER 2024)
  "TMEM155" = "SMIM43",
  "FLJ44635" = "NHSL2",
  "WI2-2373I1.2" = "FOXL3-OT1",
  "GS1-259H13.2" = "TMEM225B",
  "ZNF664-FAM101A" = "ZNF664-RFLNA",
  "MGC57346-CRHR1" = "LINC02210-CRHR1"
)

# Merge LOC genes symbols updated manually with map_genes_to_update vector.
map_genes_to_update <- c(map_genes_to_update, loc_to_update)

rm(loc_genes, not_approved_symbol, res_loc, loc_to_update)

# ADDITIONAL 3: Genes that should be removed from the Brunello library:
# SPHAR (10638): It has been removed from NCBI because is considered the 3'UTR part of the RABA4 gene.
# OCLM (10896): It has been removed from NCBI, EMBL-EBI and HGNC. There is not enougth evidence to confirm it as an independent gene.
# CRIPAK (285464): NCBI statement: This record has been withdrawn by NCBI, after discussions with CCDS collaborators. It was decided that this locus is not an independent gene.
# FAM231A (729574): NCBI statement: This record has been withdrawn by NCBI, HGNC, and EBI staff.
# ERCC6-PGBD3 (101243544): NCBI statement: This record has been withdrawn because it is no more considered as an independent readthrough gene.

# Cambiar ID de FLJ44635 A 340527

# The remaining genes that have not a suggested symbol neither were updated remain in the same way in the official pages (2025-09-23)

# Update gene symbols
gene_symbol_check <- gene_symbol_check %>%
  
  # Set previous updated rules
  mutate(
    target_gene_symbol_updated = case_when(
      Approved == TRUE | (Approved == FALSE & is.na(Suggested.Symbol)) ~ x,
      Approved == FALSE & !is.na(Suggested.Symbol) ~ Suggested.Symbol,
      TRUE ~ NA)
  ) %>%
  
  # Change specific gene symbols updated recovered with Biomart or manually.
  mutate(
    target_gene_symbol_updated = case_when(
      x %in% names(map_genes_to_update) ~ unname(map_genes_to_update[x]),
      TRUE ~ target_gene_symbol_updated
    )
  ) %>%
  dplyr::select(x, target_gene_symbol_updated) %>%
  rename(target_gene_symbol = x)
  
# To add the updated gene symbol to the brunello raw library
brunello_raw <- merge(brunello_raw, gene_symbol_check, by = "target_gene_symbol")

# Filter genes that were removed from NCBI due to not be a real gene.
genes_to_remove <- c("SPHAR", "OCLM", "CRIPAK", "FAM231A", "ERCC6-PGBD3")
  
# 5 genes and 20 sgRNA (one gene has only 2 sgRNAs) were removed from the library.
brunello_processed <- brunello_raw %>% 
  filter(!target_gene_symbol %in% genes_to_remove) %>%
  arrange(target_gene_id)

rm(ensembl, map_genes_to_update, genes_to_remove, gene_symbol_check, brunello_raw)

# ----------------------------------------------------------------- #
# ---- 3. Knonw-essentials and non-essentials preprocessing ---- ####
# ----------------------------------------------------------------- #

# This preprocessing step was performed here because the gene list are needed to generate
# processed library files for downstream analysis.
depmap_known_ess_raw <- read.csv("/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/essentials/raw/crispr_depmap_common_essentials.csv", sep = "\t")
bagel_non_ess_raw <- read.csv("/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/non_essentials/raw/NEGv1.txt", sep = "\t")

#### ----- Knonw-essentials Depmap ------ #
# Confirm if gene symbols are updated.
knonw_ess_gene_symbol_check <- checkGeneSymbols(
  depmap_known_ess_raw$GENE_SYMBOL,
  unmapped.as.na = TRUE,
  species = "human",
  expand.ambiguous = FALSE)

# One gene has outdated annotation.
knonw_ess_gene_symbol_check %>% filter(Approved == FALSE)

# Gene symbol update
depmap_known_ess_raw$GENE_SYMBOL[depmap_known_ess_raw$GENE_ID == 128061] <- "FSAF1"

# These 4 known-essential genes from Depmap are not included in the Brunello Library (Unknown reason)
lost_known_ess <- setdiff(depmap_known_ess_raw$GENE_SYMBOL, brunello_processed$target_gene_symbol_updated)

# The previous genes were removed to get a processed list of known-essential genes.
depmap_known_ess_process <- depmap_known_ess_raw %>% 
  filter(!GENE_SYMBOL %in% lost_known_ess)
  
#### ----- Known Non-essentials BAGEL2 ------ #
knonw_non_ess_gene_symbol_check <- checkGeneSymbols(
  bagel_non_ess_raw$GENE,
  unmapped.as.na = TRUE,
  species = "human",
  expand.ambiguous = FALSE)


# Known non-essentials with outdated gene symbol
knonw_non_ess_gene_symbol_check %>% 
  filter(Approved == FALSE) %>%
  # SPHAR is not a gene and PRAMEF3 has been discontinued. They had NA in suggested symbol.
  filter(!x %in% c("SPHAR", "PRAMEF3"))

# Multi symbol genes suggested were changed according to previous information.
knonw_non_ess_gene_symbol_check$Suggested.Symbol[knonw_non_ess_gene_symbol_check$x == "GIF"] <- "CBLIF"
knonw_non_ess_gene_symbol_check$Suggested.Symbol[knonw_non_ess_gene_symbol_check$x == "LOR"] <- "LORICRIN"

bagel_non_ess_process <- knonw_non_ess_gene_symbol_check %>% 
  # DO NOT exist.
  filter(!x %in% c("SPHAR", "PRAMEF3")) %>%
  dplyr::select(-Approved) %>%
  rename(GENE = x) %>%
  # Merge with BAGEL2 non-essential initial information.
  merge(bagel_non_ess_raw, ., by = "GENE") %>%
  # Keep Non-Essential genes present in the Brunello Library.
  filter(Suggested.Symbol %in% brunello_processed$target_gene_symbol_updated) %>%
  # Remove previous GENE coloum
  dplyr::select(-GENE) %>%
  # Relocate and update new GENE coloum
  relocate(Suggested.Symbol, .before = HGNC_ID) %>%
  rename(GENE = Suggested.Symbol) %>%
  # Remove duplicates
  distinct()

# Check all known-essential and known-non-essential are included in the Brunello Library
all(depmap_known_ess_process$GENE_SYMBOL %in% brunello_processed$target_gene_symbol_updated)
all(bagel_non_ess_process$GENE %in% brunello_processed$target_gene_symbol_updated)

# The final set sizes:
# DEPMAP KNOWN-ESSENTIAL GENES: 1488
# BAGEL2 NON-ESSENTIAL GENES: 911
rm(knonw_ess_gene_symbol_check, knonw_non_ess_gene_symbol_check, lost_known_ess, bagel_non_ess_raw, bagel_non_ess_raw, depmap_known_ess_raw)

# ---------------------------------------------------------------------- #
# ---- 4. Generate library files for different bioinformatics tools ---- ####
# ---------------------------------------------------------------------- #
#### ------ 3.1 MAGECK ------ ####
# Create a Brunello library file matching MAGECK requeriments.
brunello_mageck <- brunello_processed %>%
  
  # Select updated gene symbols and sgRNA sequence
  dplyr::select(target_gene_symbol_updated, sg_rna_target_sequence) %>%
  
  # Create a sgRNA ID variable with row number and gene name.
  mutate(sgRNAID = paste0(target_gene_symbol_updated, "_", row_number())) %>%
  
  # Rename variable names to match MAGECK requeriments.
  rename(Gene = target_gene_symbol_updated,
         Seq = sg_rna_target_sequence) %>%
  dplyr::select(sgRNAID, Seq, Gene) %>% 
  
  # Avoud withespaces in gene and sgRNA ID names.
  mutate(across(c(sgRNAID, Gene), ~ gsub(" ", "_", .)))


#### ------ 3.2 MAGECK OFF-TARGET sgRNAs ------ ####
# Select off-targets sgRNA IDs.
# Objective: This file is useful to normalize by off-targets sgRNA in MAGECK.
off_targets <- brunello_mageck %>%
  filter(grepl("Control", Gene)) %>%
  dplyr::select(sgRNAID)


#### ------ 3.3 MAGECK KNOWN NON-essential sgRNAs ------ ####
# Select BAGEL2 Known non-essential sgRNA IDs.
# Objective: This file is useful to normalize by known-non essentials genes in MAGECK.
known_non_essentials <- brunello_mageck %>%
  filter(Gene %in% bagel_non_ess_process$GENE) %>%
  dplyr::select(sgRNAID)

#### ------ 3.4 CRISPRcleanR ------ ####
# Objective: This Brunello library format is needed to run CRISPRcleanR.
brunello_CRISPRcleanR <- brunello_processed %>%
  
  # Create CHRM variable from NCBI Ref Code. 
  mutate(CHRM = recode(gsub("^0+", "", str_extract(genomic_sequence, "(?<=NC_0000)\\d+")), `23` = "X", `24` = "Y")) %>%
  
  # To remove "Non-Targeting Control" sgRNAs.
  filter(target_gene_symbol_updated != "Non-Targeting Control") %>%
  
  # Remove "sgRNA" from sgRNA names
  mutate(
    
    # CODE = sgRNA name + row number
    CODE = paste0(target_gene_symbol_updated, "_", row_number()),
    
    # sgRNA ENDpos is absent in the current Brunello library but it was added
    # by summing 20 to the START position.
    ENDpos = as.integer(position_of_base_after_cut_1_based + 20),
    
    # Codify strand as sense --> + and antisense --> -
    strand = case_when(strand == "sense" ~ "+",
                       strand == "antisense" ~ "-",
                       TRUE ~ NA_character_),
    
    # Add the prefix "ex" to exone information
    exon_number = paste0("ex", exon_number)
  ) %>%
  
  # Rename genes according to the CRISPRclean R requeriments.
  rename(GENES = target_gene_symbol_updated,
         EXONE = exon_number,
         STRAND = strand,
         seq = sg_rna_target_sequence,
         STARTpos = position_of_base_after_cut_1_based) %>%
  
  # Arrange by Target.Gene.id
  arrange(target_gene_id) %>%
  
  # Reorder variables in the order required by CRISPRcleanR
  dplyr::select(CODE, GENES, EXONE, CHRM, STRAND, STARTpos, ENDpos, seq) 


#### ------ 3.5 SSC ------ ####
# Objective: This Brunello library format is needed to run SSC tool to generate
# sgRNA effiencies to use in MAGECK-MLE.
brunello_SSC <- brunello_mageck %>%
  relocate(Seq, .before = sgRNAID) %>%
  rename(Sequence = Seq)


head(brunello_CRISPRcleanR)
head(brunello_mageck)
head(depmap_known_ess_process)
head(bagel_non_ess_process)
head(known_non_essentials)
head(off_targets)
head(brunello_SSC)

# PRIMERA MRD.

#### ------ 3.6 BAGEL2 PREPARE ALIGNMENT FILE------ ####
# Objective: BAGEL2 has an optional argument to use an aligment file which provide
# information about the mapping of each sgRNA with additional protein-coding genes,
# non-coding regions... This file is obtained after align sgRNA sequences with the
# human genome using bowtie1 and a custom function created by BAGEL2 developers.
# This function requires a library file following the next structure:
# READID: sgRNA ID
# GUIDESEQ: sgRNA sequence
bagel2_align_lib <- brunello_mageck %>%
  dplyr::select(-Gene) %>%
  rename(READID = sgRNAID,
         GUIDESEQ = Seq)

# Here we add other Brunello formats needed...


#### ------ 4. SAVE PROCESSED BRUNELLO LIBRARIES ----- ####
# write.table(brunello_CRISPRcleanR, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/brunello/processed/v2/brunello_library_CRISPRcleanR.txt", sep = "\t", quote = F, row.names = F)
# write.table(brunello_mageck, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/brunello/processed/v2/brunello_library_MAGECK.txt", sep = "\t", quote = F, row.names = F)
# write.table(off_targets, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/brunello/processed/v2/brunello_library_off_target_sgRNAs.txt", sep = "\t", quote = F, row.names = F)
# write.table(known_non_essentials, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/brunello/processed/v2/brunello_library_non_essentials_sgRNAs.txt", sep = "\t", quote = F, row.names = F)

#write.table(bagel2_align_lib, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/brunello/processed/v2/brunello_library_bagel2_alignment.txt", sep = "\t", quote = F, row.names = F)

# Additional V2
# write.table(brunello_SSC, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/sgRNA_libraries/processed/v1/brunello_library_SSC.in", sep = "\t", quote = F, row.names = F)

#### ------ 4. SAVE GENE LISTS ----- ####
# write.table(depmap_known_ess_process, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/essentials/processed/crispr_depmap_common_essentials_processed.tsv", sep = "\t", quote = F, row.names = F)
# write.table(bagel_non_ess_process, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/non_essentials/processed/common_non_essentials_processed.tsv", sep = "\t", quote = F, row.names = F)



# ---- PREPARE FUNCTIONAL AND ESSENTIAL GENE SETS ---- #
# Se evaluó el conjunto de genes considerado como HISTONES EN CRISPRcleanR pero solo un 14% eran common-esentials
# en depamp y no se incluyeron aqui.
# 1. Download KEGG LEGACY pathways from MSIGDB
msig <- msigdbr(db_species = "HS", species = "human", collection =  "C2", subcollection = "CP:KEGG_LEGACY") 


# 2. Select pathways to evaluate.
pathways <- c(
  "KEGG_DNA_REPLICATION",
  "KEGG_RNA_POLYMERASE",
  "KEGG_SPLICEOSOME",
  "KEGG_PROTEASOME",
  "KEGG_RIBOSOME",
  "KEGG_CELL_CYCLE",
  "KEGG_AMINOACYL_TRNA_BIOSYNTHESIS",
  "KEGG_PURINE_METABOLISM",
  "KEGG_PYRIMIDINE_METABOLISM"
)


# 3. Extract genes per pathway of interest.
pathway_list <- setNames(
  lapply(pathways, function(p) {
    msig %>% filter(gs_name == p) %>% pull(gene_symbol) %>% unique()
  }),
  pathways
)

# 4. Check outdated gene symbols in KEGG pathway genes.
pathway_checks <- lapply(pathway_list, function(p) {
  checkGeneSymbols(p, unmapped.as.na = TRUE, species = "human", expand.ambiguous = FALSE)}
  )

# All symbols are updated in all KEGG pathways.
lapply(pathway_checks, function(p) {
  table(p$Approved, useNA = "always")}
)

# 5. To extract common essentials depmap in each pathway
pathways_depmap <- lapply(pathway_list, function(p) {
  intersect(p, depmap_known_ess_process$GENE_SYMBOL)
})


# 6. Assesment.

# 6.1 Initial pathway size
lapply(pathway_list, function(p){
  length(p)
}
)

# 6.2 Intersect pathway and Depmap Core Essentials
lapply(pathway_list, function(p){
  length(intersect(p, depmap_known_ess_process$GENE_SYMBOL)) / length(p)
}
)

# 6.3 Number of non-known essential depmap genes in the pathway
lapply(pathway_list, function(p){
  length(setdiff(p, depmap_known_ess_process$GENE_SYMBOL))
}
)

# Selected list for further use
pathways_selected <- c(
  "KEGG_DNA_REPLICATION",
  "KEGG_RNA_POLYMERASE",
  "KEGG_SPLICEOSOME",
  "KEGG_PROTEASOME",
  "KEGG_RIBOSOME"
)


# 6. Final lists
# Raw pathways
df_pathways_selected <- lapply(pathways_depmap, function(genes) {data.frame(GENE_SYMBOL = genes)})
df_pathways_selected <- df_pathways_selected[pathways_selected]

# Save the raw and processed data.
# Raw
save(pathway_list, file = "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/functional_sets/raw/essential_functional_sets.RData")
pathway_list_raw <- lapply(pathway_list, function(genes) {data.frame(GENE_SYMBOL = genes)})
writexl::write_xlsx(pathway_list_raw, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/functional_sets/raw/essential_functional_sets.xlsx")

# Processed
save(df_pathways_selected, file = "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/functional_sets/processed/essential_functional_sets_processed.RData")
writexl::write_xlsx(df_pathways_selected, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/gene_lists/functional_sets/processed/essential_functional_sets_processed.xlsx")



# ---------------------------------------------------------------------- #
# ---------------- 5. STRING V12 NETWORK FILE ----------------------- ####
# ---------------------------------------------------------------------- #
# Full STRING interaction file and protein info file.
string_network_raw <- read.delim("/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/network_files/STRING/raw/9606.protein.links.full.v12.0.txt", sep = " ")
string_prot_info_raw <- read.delim("/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/network_files/STRING/raw/9606.protein.info.v12.0.txt", sep = "\t")

# Pre-process
colnames(string_prot_info_raw) <- gsub("X.", "", colnames(string_prot_info_raw))
string_prot_info_raw <- string_prot_info_raw %>% dplyr::select(string_protein_id, preferred_name)


string_network_processed <- string_network_raw %>% 
  
  # Filter high-confidence interactions --> >700 according to STRING information.
  filter(combined_score >= 700) %>%
  
  # Select protein1 and protein2
  dplyr::select(protein1, protein2) %>%
  
  # Map proteins in "protein1" wich their gene symbol.
  left_join(string_prot_info_raw, by = c("protein1" = "string_protein_id")) %>%
  rename(protein1_symbol = preferred_name) %>%
  
  # Map proteins in "protein2" wich their gene symbol.
  left_join(string_prot_info_raw, by = c("protein2" = "string_protein_id")) %>% 
  rename(protein2_symbol = preferred_name) %>% 
  
  # Extract gene symbols pairs.
  dplyr::select(protein1_symbol, protein2_symbol) %>%
  
  # Detect duplicate interactions due to reverse order in the pair_canonical variable
  mutate(pair_canonical = map2_chr(protein1_symbol, protein2_symbol, ~paste(sort(c(.x, .y)), collapse = "_")))
  

# Check duplicated interactions. All interactions are duplicated, so we only need one pair.
duplicated_interactions <- string_network_processed %>% filter(duplicated(pair_canonical) | duplicated(pair_canonical, fromLast = T)) %>% arrange(pair_canonical)

length(unique(string_network_processed$pair_canonical))*2 == nrow(string_network_processed)


string_network_processed <- string_network_processed %>% 
  
  # Duplicates removal
  distinct(pair_canonical, .keep_all = TRUE) %>% 
  dplyr::select(protein1_symbol, protein2_symbol) %>%

  # Remove interactions which cotains ENSEMBL Protein or Uniprot codes because they are not in Brunello Library.
  filter(!(
    str_detect(protein1_symbol, "^ENSP") |
    str_detect(protein1_symbol, "_HUMAN$") |
    str_detect(protein2_symbol, "^ENSP") |
    str_detect(protein2_symbol, "_HUMAN$")
  
  ))


setdiff(unique(string_network_processed$protein1_symbol), unique(brunello_processed$target_gene_symbol_updated))

# Protein1 CHECKING
# Check gene symbol to know if they are updated
protein1_symbol_check <- checkGeneSymbols(
  unique(string_network_processed$protein1_symbol),
  unmapped.as.na = TRUE,
  species = "human",
  expand.ambiguous = FALSE)


protein1_symbol_check <- protein1_symbol_check %>%
  
  # Set previous updated rules
  mutate(
    protein1_symbol_updated = case_when(
      Approved == TRUE | (Approved == FALSE & is.na(Suggested.Symbol)) ~ x,
      Approved == FALSE & !is.na(Suggested.Symbol) ~ Suggested.Symbol,
      TRUE ~ NA)
  ) %>%
  dplyr::select(x, protein1_symbol_updated) %>%
  rename(protein1_symbol = x)


# Update cases where double symbol is found according to Brunello Library.
protein1_symbol_check <- protein1_symbol_check %>%
  mutate(
    protein1_symbol_updated = case_when(
      protein1_symbol == "MPP6" ~ "PALS2",
      protein1_symbol == "LOR"  ~ "LORICRIN",
      protein1_symbol == "TAZ"  ~ "WWTR1",
      TRUE ~ protein1_symbol_updated
    )
  )


# Protein2 CHECKING
# Check gene symbol to know if they are updated
protein2_symbol_check <- checkGeneSymbols(
  unique(string_network_processed$protein2_symbol),
  unmapped.as.na = TRUE,
  species = "human",
  expand.ambiguous = FALSE)


# PROTEIN2 CHECK SYMBOL
protein2_symbol_check <- protein2_symbol_check %>%
  
  # Set previous updated rules
  mutate(
    protein2_symbol_updated = case_when(
      Approved == TRUE | (Approved == FALSE & is.na(Suggested.Symbol)) ~ x,
      Approved == FALSE & !is.na(Suggested.Symbol) ~ Suggested.Symbol,
      TRUE ~ NA)
  ) %>% 
  dplyr::select(x, protein2_symbol_updated) %>%
  rename(protein2_symbol = x)


# Update cases where double symbol is found according to Brunello Library.
protein2_symbol_check <- protein2_symbol_check %>%
  mutate(
    protein2_symbol_updated = case_when(
      protein2_symbol == "MPP6" ~ "PALS2",
      protein2_symbol == "LOR"  ~ "LORICRIN",
      protein2_symbol == "TAZ"  ~ "WWTR1",
      TRUE ~ protein2_symbol_updated
    )
  )



# Update gene symbols in the network file

string_network_processed <- string_network_processed %>% 
  left_join(protein1_symbol_check, by = c("protein1_symbol" = "protein1_symbol")) %>% 
  left_join(protein2_symbol_check, by = c("protein2_symbol" = "protein2_symbol")) %>% 
  dplyr::select(protein1_symbol_updated, protein2_symbol_updated)


# Some genes in the network file are not included in the Brunello Librarl.
setdiff(unique(string_network_processed$protein2_symbol_updated), unique(brunello_processed$target_gene_symbol_updated))%>% sort()
setdiff(unique(string_network_processed$protein1_symbol_updated), unique(brunello_processed$target_gene_symbol_updated)) %>% sort()

#' NOTE (2025-10-15): Genes not included in the Brunello library were not excluded from the network file in this version because
#' it is not clear how BAGEL2 uses this network file. The question is if BAGEL2 uses all interactions of a gene without notify if the paired
#' is included or not in the CRISPR screening library.
#' Example: prot1 (in Brunello) - prot2 (not in Brunello) --> is it useful?
length(unique(string_network_processed$protein2_symbol_updated))
length(unique(brunello_processed$target_gene_symbol_updated))

rm(string_prot_info_raw, string_network_raw, protein1_symbol_check, protein2_symbol_check, duplicated_interactions)


# SAVE processed file.
#write.table(string_network_processed, "/home/andres/Escritorio/bioinfo_hematology_unit/crispr_hunt/workflow/resources/network_files/STRING/processed/human_protein_network_STRING_processed.txt", sep = "\t", quote = F, row.names = F)


   