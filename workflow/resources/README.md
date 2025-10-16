# CRISPR-Hunt Resources

Authoritative reference datasets for gene essentiality, non-essentiality, and sgRNA efficiency. These resources are versioned, documented, and designed for direct use in CRISPR screening pipelines.

## sgRNA libraries

### Brunello - raw
- **brunello_library_raw_v1.txt**:
  - Source: Addgene
  - Link: https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt
  - Description: Brunello Library Target Genes. It was a version 1 because the decimal separator in the variable "Rule Set 2 score" was not properly recorded.
  - Download date: 2025-07-01
 
- **brunello_library_raw_v2.txt**:
  - Source: Addgene
  - Link: https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt
  - Description: Brunello Library Target Genes. It was the corrected version from the v1 file.
  - Download date: 2025-09-24

### Brunello - processed
- Version 1 (v1)
   - Raw file: `broadgpp_brunello_library_raw_v1.txt`
   - Files included: 
   	- `brunello_library_CRISPRcleanR.txt`: Processed Brunello library formatted for use with CRISPRcleanR (N = 76 442).
   	- `brunello_library_MAGECK.txt`: Processed Brunello library formatted for use with MAGECK (N = 77 442).
   	- `brunello_library_non_essentials_sgRNAs.txt`: sgRNA IDs for non-essential genes (BAGEL2), generated from the NEGv1.txt raw file. Used for MAGECK normalization (N = 3 584).
   	- `brunello_library_off_target_sgRNAs.txt`: sgRNA IDs for off-target sgRNAs. Used for MAGECK normalization (N = 1 000).
   	- `brunello_library_bagel2_alignment.txt`: Processed Brunello library formatted for use with precalc_library_alignment_info.py as a prior BAGEL2 step to identify multi-target sgRNA effects.
  - Generated: 2025-07-23 (SSC-formatted file created on 2025-09-24).
  - Preprocessing script: These files were created in the Linux Shell.
  - Description: Preprocessed Brunello Library providing multiple formats and supporting files for MAGECK and CRISPRcleanR analyses. Supporting files include sgRNA IDs required for normalization by off-target sgRNAs and non-essential genes in MAGECK.
	
- Version 2 (v2)
  - Raw file: `broadgpp_brunello_library_raw_v2.txt`
  - Files included: Same as v1. 
  - Generated: 2025-09-24
  - Preprocessing script: Brunello_library_preprocessing.R
  - Changes: Updated gene symbols; obsolete or non-existent genes were removed. After preprocessing, the total number of sgRNAs decreased to 77,423. The number of known non-essential sgRNAs increased to 3,650.

- Additional file generated from raw v2:
  - brunello_library_SSC.in (Generated 2025-09-24) — First-time generation of SSC-formatted Brunello library from raw v2.


## FastQC
- **adapter_list.txt**: 
  - Source: FastQC Github.
  - Link: https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt
  - Description: Adapter file used in FastQC to detect primers in input sequences. 5'- and 3'- adapter from Brunello CRISPR screening library were added into the list.
  - Download date: 2025-06-01

## Gene Lists

### Essential genes - Raw
- **CEGv2.txt**
  - Source: BAGEL2 GitHub
  - Link: https://github.com/hart-lab/bagel/blob/master/CEGv2.txt
  - Version: v2 (Second version in BAGEL2 Github)
  - Description: Core essential genes
  - Download Date: 2025-08-08

- **pan_species_common_essentials.txt**:
  - Source: BAGEL2 GitHub 
  - Link: https://github.com/hart-lab/bagel/blob/master/pan-species-control-essentials-50genes.txt
  - Version: v1 (First version in BAGEL2 Github)
  - Original name: pan-species-control-essentials-50genes.txt
  - Description: Core essential genes found in Human, Rat and Mouse.
  - Download Date: 2025-08-08

- **crispr_depmap_common_essentials.csv**
  - Source: DepMap 
  - Link: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=CRISPRInferredCommonEssentials.csv
  - Version: DepMap, Broad (2025). DepMap Public 25Q2. 
  - Original name: CRISPRInferredCommonEssentials.csv
  - Download Date: 2025-08-08 
  - Description: List of genes identified as dependencies across all lines (Post-Chronos). Total: 1492

- **achilles_common_essentials.csv**
  - Source: DepMap
  - Link: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=AchillesCommonEssentialControls.csv
  - Version: DepMap, Broad (2025). DepMap Public 25Q2. 
  - Original name: AchillesCommonEssentialControls.csv
  - Download Date: 2025-08-08
  - Description: List of genes used as positive controls, intersection of Biomen (2014) and Hart (2015) essentials (Pre-Chronos).


### Essential genes - Processed
- **common_essentials_processed.tsv**
  - Raw file: `crispr_depmap_common_essentials.csv`
  - Processing: Litle changes. One gene symbol was updated and 4 genes were removed to not be in the Brunello Library.
  - Generated: 2025-09-26 
  - Preprocessing script: Brunello_library_preprocessing.R
  - Description: Depmap (Post-Chronos) Core essential genes identified across all lines processed for CRISPR screening analysis. Total: 1488 genes.

### Non-Essential genes - Raw
- **NEGv1.txt**
  - Source: BAGEL2 GitHub
  - Link: https://github.com/hart-lab/bagel/blob/master/NEGv1.txt
  - Version: v1 (First version in BAGEL2 Github)
  - Description: Core non-essential genes (Hart 2015). Total: 927 genes.
  - Download Date: 2025-08-08

- **pan_species_non_essentials.txt**:
  - Source: BAGEL2 GitHub
  - Link: https://github.com/hart-lab/bagel/blob/master/pan-species-control-nonessentials-50genes.txt
  - Version: v1 (First version in BAGEL2 Github)
  - Original name: pan-species-control-nonessentials-50genes.txt
  - Description: Core non-essential genes found in Human, Rat and Mouse.
  - Download Date: 2025-08-08

- **achilles_non_essentials.csv**
  - Source: DepMap
  - Link: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=AchillesNonessentialControls.csv
  - Version: DepMap, Broad (2025). DepMap Public 25Q2.
  - Original name: AchillesNonessentialControls.csv
  - Download Date: 2025-08-08
  - Description: List of genes used as negative controls (Hart (2014) nonessentials) (Pre-Chronos).

### Non-Essential genes - Processed
- **common_non_essentials_processed.tsv**
  - Raw file: `NEGv1.txt`
  - Processing: Gene symbols were updated and genes not preseNt in the Brunello Library were removed.
  - Generated: 2025-09-26
  - Preprocessing script: Brunello_library_preprocessing.R
  - Description: Processed core non-essential genes (Hart 2015). Total 911 genes.

### Functional sets - Raw
- **essential_functional_sets**:
  - Source: MSigDB 2025.1
  - Link: https://www.gsea-msigdb.org/gsea/index.jsp
  - Description: Complete initial download of pathways considered as potential "essential pathways," based on literature and the CRISPRcleanR paper, from KEGG. The pathways were retrieved using the msigdbr (v25.1.1) R package. Included pathways:
  	- KEGG_DNA_REPLICATION
  	- KEGG_RNA_POLYMERASE
  	- KEGG_SPLICEOSOME
  	- KEGG_PROTEASOME
  	- KEGG_RIBOSOME
  	- KEGG_CELL_CYCLE
  	- KEGG_AMINOACYL_TRNA_BIOSYNTHESIS
  	- KEGG_PURINE_METABOLISM
  	- KEGG_PYRIMIDINE_METABOLISM.
  - Script: Brunello_library_preprocessing.R (⚠️ TODO: migrate to an independent script)
  - Formats: .RData and .xlsx.
  - Download Date: 2025-10-10

### Functional sets - Processed
- **essential_functional_sets_processed**
  - Raw file: `essential_functional_sets_raw.RData`
  - Processing: Each pathway was filtered to retain genes overlapping with DepMap known-essential genes. Only pathways with >65% overlap were kept, and only the corresponding DepMap known-essential genes for each pathway were included. Final pathways:
    	- KEGG_DNA_REPLICATION
  	- KEGG_RNA_POLYMERASE
  	- KEGG_SPLICEOSOME
  	- KEGG_PROTEASOME
  	- KEGG_RIBOSOME
  - Generated: 2025-10-10
  - Script: Brunello_library_preprocessing.R (⚠️ TODO: migrate to an independent script)
  - Formats: .RData and .xlsx.
  - Description: Curated functional gene sets suitable for evaluating CRISPR screening quality control (QC).

## sgRNA Efficiency Files

This directory contains files related to sgRNA efficiency for the Brunello library, used in MAGECK MLE calculations.

- **brunello_library_sgRNA_efficiency_rule_set_2.txt**
  - Source: Addgene Brunello Library
  - Metric: Rule Set 2 score
  - Description: This metric is an value for sgRNA efficiency created by Brunello developers .Initially tested as input for MAGECK MLE. It is not the recommended sgRNA efficiency information by MAGECK developers. It was not used in final MAGECK MLE analysis but it could be used to map sgRNA effiency in every experiment.
  - Download Date: 2025-07-01

- **brunello_library_sgRNA_efficiency_SSC.out**
  - Generated from: Preprocessed Brunello library (`brunello_library_SSC.in`)
  - Tool: SSC (Source: https://sourceforge.net/projects/spacerscoringcrispr/)
  - Description: sgRNA efficiency predicted following MAGECK tutorial recommendations (SCC tool). This file was used as input for MAGECK MLE.
  - Generation Date: 2025-09-24 
  - Generation steps:
    1. SSC tool downloaded (https://sourceforge.net/projects/spacerscoringcrispr/files/SSC0.1/) and compiled using `make` in the project directory and README file instructions.
    2. Input library preprocessed with `Brunello_library_preprocessing.R`. The initial Brunello Library was processed to create the required SSC format.
    3. SSC executed with command:
       ```
       ./SSC -l 20 -m SSC0.1/matrix/human_CRISPRi_20bp.matrix \
            -i <input_library> \
            -o <output_file>
       ```
       - `-l`: input sequence length (20 bp) because is the sgRNA sequence lenght in the Brunello Library.
       - `-m`: matrix matching sequence length (human_CRISPRi_20bp.matrix). This file is recommended in the MAGECK MLE advanced tutorial.
       - `-i`: preprocessed Brunello library.
       - `-o`: output efficiency file formatted for MAGECK MLE
 
### genomes

- **human_GRC38/GRCh38.primary_assembly.genome.fa.gz**.
  - Source: EMBL-EBI
  - Link: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
  - Version: GENCODE Human Release 49 (2025-09-02)
  - Description: This FASTA file corresponds to the GRCh38 primary assembly provided by the GENCODE consortium. It was selected as the reference genome for the generation of the sgRNA alignment information used in BAGEL2 multi-target filtering.
The Brunello CRISPR library (Doench et al., 2016) was originally designed around 2015, and no explicit information about the genome build or annotation used in its design is provided in the original publication and no information was found in other web sites. Manual alignment of Brunello sgRNAs against current human genome assemblies confirmed that their genomic coordinates are consistent with GRCh38 (hg38).
Given that only minor sequence corrections have been introduced in subsequent GRCh38 patch releases (mostly affecting alternative loci and non-coding regions), the most recent GRCh38 primary assembly was adopted (GRCh38.p14). This choice ensures:

		- Compatibility with modern CRISPR analysis tools (e.g., BAGEL2, MAGeCK, CRISPick).

		- Full consistency with the paired GENCODE Release 49 annotation (GTF).

		- Reproducibility and forward-compatibility for future analyses.
In summary, this genome build provides an accurate and up-to-date reference for sgRNA alignment while maintaining backward consistency with the coordinates of the Brunello library.
  - Download Date: 2025-10-14
  
- **human_GRC38/gencode.v49.primary_assembly.annotation.gtf**:
  - Source: EMBL-EBI
  - Link: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
  - Version: GENCODE Human Release 49 (2025-09-02)
  - Description: GFT annotation file for the GRCh38.primary_assembly.genome.fa.gz FASTA file.
  - Download Date: 2025-10-14

- **human_GRC38/GRCh38_primary_bt1_index**
  - Generated from: GRCh38.primary_assembly.genome.fa.gz
  - Indexing Tool: Bowtie v1.3.1 (downloaded inside bagel2 mamba environment).
  - Command line: bowtie-build GRCh38.primary_assembly.genome.fa GRCh38_primary_bt1_index.
  - Description: This folder contains the Bowtie1 genome index built from the GENCODE GRCh38 primary assembly (release 49). The index is required by the precalc_library_alignment_info.py script to align each sgRNA sequence to the reference genome and generate the alignment-info file used by BAGEL2 for multi-target guide filtering. Bowtie version 1 (not Bowtie2) was explicitly selected because the script and downstream analysis expect the alignment behavior of Bowtie1, which is optimized for short, exact reads (~20 bp), as used in CRISPR screening libraries.
  - Contents:

	- GRCh38_primary_bt1_index.1.ebwt

	- GRCh38_primary_bt1_index.2.ebwt

	- GRCh38_primary_bt1_index.3.ebwt

	- GRCh38_primary_bt1_index.4.ebwt

	- GRCh38_primary_bt1_index.rev.1.ebwt

	- GRCh38_primary_bt1_index.rev.2.ebwt
  - Log file: GRCh38_primary_bt1_index.log.
  - Date of generation: 2025-10-14

### alignment_info

- **alignment_info/brunello_bagel2_alignment_info.txt**
  - Generated from:
	- Genome index: human_GRCh38/GRCh38_primary_bt1_index (Bowtie v1 format)
	- Annotation: human_GRCh38/gencode.v49.primary_assembly.annotation.gtf
	- Input library: sgRNA_libraries/brunello/processed/v2/brunello_library_bagel2_alignment.txt
  - Tools:
  	- Bowtie v1.3.1 (using in the py script)
	- Python script: precalc_library_alignment_info.py (BAGEL2 pre-processing utility)
  - Description: This file contains precomputed genome alignment information for the Brunello CRISPR library. Each sgRNA sequence is aligned against the human reference genome (GRCh38) to identify potential off-target matches and multi-target effects. The resulting table is used by BAGEL2 to exclude sgRNAs with non-specific or ambiguous targeting efficiency.
  - Command line (executed within the bagel2 mamba environment): 
  	```python3 path/to/precalc_library_alignment_info.py 
  	path/to/brunello_library_bagel2_alignment.txt 
  	path/to/human_GRCh38/GRCh38_primary_bt1_index \
  	path/to/GRCh38_primary_bt1_index human_GRCh38/gencode.v49.primary_assembly.annotation.gtf 
  	--output path/to/alignment_info/brunello_bagel2_alignment_info.txt 
  	--custompam NGG 
  	--pam-loc 0 
  	2> path/to/alignment_info/brunello_bagel2_alignment_info.log```
  - Log file: brunello_bagel2_alignment_info.log captures detailed Bowtie and parsing diagnostics.
  - Generation Date: 2025-09-24 
  
  
### network_files - raw

- **STRING/raw/9606.protein.links.full.v12.0.txt**
  - Source: STRING Database
  - Link: https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz
  - Version: v12.0
  - Description: Comprehensive human protein–protein interaction network from STRING. This file includes all interactions and each interaction is scored both per evidence channel and as a combined score, reflecting the confidence of the interaction.
  - Download Date: 2025-10-15

- **STRING/raw/9606.protein.info.v12.0.txt**
  - Source: STRING Database
  - Link: https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
  - Version: v12.0
  - Description: List of STRING proteins including their display names and descriptions.
  - Download Date: 2025-10-15
  
### network_files - processed

- **STRING/processed/human_protein_network_string_processed.txt**
  - Generated from:
	- STRING full protein links: 9606.protein.links.full.v12.0.txt (STRING human protein–protein interaction network)
	- STRING protein : 9606.protein.info.v12.0.txt
  - Script: Brunello_library_preprocessing.R (⚠️ TODO: migrate to an independent script)
  - Processing: 
  	- High-confidence filtering: retained only interactions with combined_score >= 700, corresponding to high confidence according to STRING scoring criteria.
	- Symbol mapping: converted STRING protein IDs (Ensembl Protein IDs) to gene symbols using the mapping provided in protein.info.
	- Duplicate removal: interactions appearing in both orientations (e.g. A–B and B–A) were deduplicated to retain a single representative entry.
	- Non-matching entries cleanup: removed interactions in which at least one interactor was identified only by an Ensembl or UniProt code (i.e. no corresponding gene symbol).
	- Gene symbol harmonization: updated all gene symbols to match Brunello CRISPR library nomenclature, following the same mapping and curation logic applied during library preprocessing.
  - Description: High-confidence, curated human protein–protein interaction network compatible with BAGEL2 analysis. Number of protein interactions: 235,606
  - Generation Date: 2025-10-15
  
