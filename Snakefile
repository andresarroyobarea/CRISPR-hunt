import glob
import os
import sys
import pandas as pd
from workflow.utils.common import (
    get_samples,
    get_resource,
    get_count_table_to_norm,
    get_count_table_to_test,
    get_lfc_file,
)
from snakemake.utils import min_version

# ---- Snakemake minimal version ---- #
min_version(
    "9.1.1"
)  # Once the pipeline is complete, previous Snakemake versions should be checked.


# ---- Read Config file ---- #
configfile: "config/config.yaml"


# ---- Set up report ---- #
# TODO: This section will cover the generation of Snakemake reports.


# ---- Docker image ---- #
# TODO: A container to define the underlying OS for each job when using the Snakemake workflow with --use-conda --use-singularity.
# Example: container: "docker://continuumio/miniconda3"

# ---- Config and global variables ---- #
# Project name
project = config["project"]

# Samples
samples = config["samples"]
sample_filters = config["sample_filters"]

df_samples = pd.read_csv(samples, sep="\t")

SAMPLES = get_samples(
    df_samples,
    included=sample_filters.get("included", True),
    batch=sample_filters.get("batch"),
    condition=sample_filters.get("condition"),
)

# print("Muestras seleccionadas:",SAMPLES)

# Conditions
# Filtered samples
df_samples_fil = df_samples[df_samples["sample"].isin(SAMPLES)]

# treated samples
treat_cond = config["parameters"]["mageck_rra_test"]["treat_cond"]

# control samples
ctrl_cond = config["parameters"]["mageck_rra_test"]["ctrl_cond"]

treat_samples = (
    df_samples_fil[df_samples_fil["condition"] == treat_cond]
    .sort_values("replicate")["sample"]
    .tolist()
)
ctrl_samples = (
    df_samples_fil[df_samples_fil["condition"] == ctrl_cond]
    .sort_values("replicate")["sample"]
    .tolist()
)

# ---- Design matrix ---- #
# To select filtered samples from the original desing matrix.
design = config["design"]
design_matrix = pd.read_csv(design, sep="\t", index_col=False)

# print(design_matrix.columns)

design_matrix_filt = design_matrix[design_matrix["sample"].isin(SAMPLES)]
design_matrix_filt_path = "workflow/resources/design_matrix/design_matrix_filtered.tsv"
design_matrix_filt.to_csv(design_matrix_filt_path, sep="\t", index=False)

# print(df_samples_fil)
# print(design_matrix_filt)

# Normalization status for BAGEL2
# bagel2_norms = ensure_raw_status(config.get("bagel2_norms", [""]))
bagel2_norms = config.get("bagel2_norms", [""])

# Normalization status for MAGECK count
mageck_norms = config.get("mageck_norms", [""])

# print(mageck_norms)

# MAGECK RRA / MLE tests
mageck_test = [n for n in mageck_norms if n != "raw"]

if config.get("cnv_correction", {}).get("enabled", False):
    mageck_test.append("cnvcorr")

# CNV impact comparisons
# uncorr_norm = config["parameters"]["cnv_impact"][""]


# ---- RULE MODULES ---- #
include: "workflow/rules/qc.smk"


if config["parameters"]["trimming"].get("enabled", False):

    include: "workflow/rules/trimming.smk"


include: "workflow/rules/mageck_count_raw.smk"
include: "workflow/rules/filter_sgrna_counts.smk"
include: "workflow/rules/mageck_normalize.smk"


if config["cnv_correction"].get("enabled", False):

    include: "workflow/rules/cnv_correction.smk"
    # include: "workflow/rules/cnv_correction_qc.smk"
    include: "workflow/rules/cnv_impact.smk"


for method in config.get("essentiality_methods", []):
    if method == "mageck_rra":

        include: "workflow/rules/mageck_rra_test.smk"

    elif method == "mageck_mle":

        include: "workflow/rules/mageck_mle_test.smk"

    elif method == "bagel2":

        include: "workflow/rules/bagel2.smk"


# ---- Target rules ---- #
# TODO: This is a initial version to collect the final MAGECK RRA output but it will
# be updated for general purposes.
rule all:
    input:
        "results/qc_raw/multiqc_report.html",
        "results/qc_trimming/multiqc_report.html",
        expand(
            "results/qc_trimming/{sample}/{sample}_trimmed_fastqc.html", sample=SAMPLES
        ),
        expand(
            "results/mageck_count_raw/{project}.count.txt",
            project=project,
        ),
        expand(
            "results/filter_sgrna_counts/{project}_processed.count.txt",
            project=project,
        ),
        expand(
            "results/mageck_normalize/{mageck_norms}/{project}_{mageck_norms}.count_normalized.txt",
            project=project,
            mageck_norms=mageck_norms,
        ),
        expand(
            "results/mageck_rra_NCT_qc/{project}_processed_NCT.sgrna_summary.txt",
            project=project,
        ),
        expand(
            "results/mageck_rra_NCT_qc/{project}_processed_NCT.gene_summary.txt",
            project=project,
        ),
        expand(
            "results/mageck_rra_test/{mageck_test}/{project}_{mageck_test}.sgrna_summary.txt",
            project=project,
            mageck_test=mageck_test,
        ),
        expand(
            "results/mageck_rra_test/{mageck_test}/{project}_{mageck_test}.gene_summary.txt",
            project=project,
            mageck_test=mageck_test,
        ),
        expand(
            "results/mageck_mle_test/{mageck_test}/{project}_{mageck_test}.sgrna_summary.txt",
            project=project,
            mageck_test=mageck_test,
        ),
        expand(
            "results/mageck_mle_test/{mageck_test}/{project}_{mageck_test}.gene_summary.txt",
            project=project,
            mageck_test=mageck_test,
        ),
        #expand(
        #    "results/cnv_correction_qc_rra/{project}_cnvcorr_rra_QC.pdf",
        #    project=project,
        #),
        #expand(
        #    "results/cnv_correction_qc_mle/{project}_cnvcorr_mle_QC.pdf",
        #    project=project,
        #),
        expand(
            "results/cnv_impact/{mageck_norms}/{project}_{mageck_norms}_cnv_impact_stats.xlsx",
            project=project,
            mageck_norms=mageck_norms,
        ),
        expand(
            "results/cnv_impact/{mageck_norms}/{project}_{mageck_norms}_cnv_impact_plots.pdf",
            project=project,
            mageck_norms=mageck_norms,
        ),
        expand(
            "results/bagel2_fc/{project}_total.foldchange",
            project=project,
        ),
        expand(
            "results/bagel2_fc/{project}_total.normed_readcount",
            project=project,
        ),
        expand(
            "results/bagel2_bf/gene_level/{project}_{norm_state}_gene.bf",
            project=project,
            norm_state=bagel2_norms,
        ),
        expand(
            "results/bagel2_bf/sgrna_level/{project}_{norm_state}_sgrna.bf",
            project=project,
            norm_state=bagel2_norms,
        ),
        expand(
            "results/bagel2_pr/gene_level/{project}_{norm_state}-pr",
            project=project,
            norm_state=bagel2_norms,
        ),
        expand(
            "results/filter_sgrna_counts/{project}_processed_NCT.count.txt",
            project=project,
        ),
        expand(
            "results/cnv_correction/{project}_cnvcorr.count.txt",
            project=project,
        ),
        expand(
            "results/cnv_correction/{project}_cnvcorr.foldchange",
            project=project,
        ),
