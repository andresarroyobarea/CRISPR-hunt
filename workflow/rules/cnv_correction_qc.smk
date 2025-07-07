rule cnv_correction_qc_rra:
    input: 
        sgrna_summary="results/mageck_mle_test/{project}_sgrna_summary.txt",
        gene_summary="results/mageck_mle_test/{project}_gene_summary.txt",
    output: 
        qc_report="results/cnv_correction_qc_rra/{project}_QC.pdf",
    conda:
        config["conda_envs"]["crisprcleanR"],
    params:
        library_type=config["parameters"]["cnv_correction"]["library_type"],
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"],
        selection_type=confif["parameters"]["cnv_correction_qc"]["selection_type"],
        essential_genes=config["parameters"]["bagel2_fc"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_fc"]["non_essential_genes"],
    log:
        "logs/cnv_correction_qc/{project}_cnv_correction_qc.log",
    benchmark:
        "benchmarks/cnv_correction_qc/{project}_cnv_correction_qc.bmk",
    shell:
        """
        Rscript scripts/run_crisprcleanR_qc_rra.R
            --sgrna-summ {input.sgrna_summary} \
            --gene-summ {input.gene_summary} \
            --lib-type {params.library_type} \
            --sgrna-library {params.sgrna_library} \
            --sel-type {params.selection_type} \
            --essen-genes {params.essential_genes} \
            --non-essen-genes {params.non_essential_genes}
        
        """


rule cnv_correction_qc_mle:
    input: 
        gene_summary="results/mageck_mle_test/{project}_gene_summary.txt",
    output:
        ""
    conda:
        config["conda_envs"]["crisprcleanR"],
    params:
        library_type=config["parameters"]["cnv_correction"]["library_type"],
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"],
        essential_genes=config["parameters"]["bagel2_fc"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_fc"]["non_essential_genes"],
        label=config["project"]
    log:
        "logs/cnv_correction_qc_mle/{project}_cnv_correction_qc_mle.log",
    benchmark:
        "benchmarks/cnv_correction_qc_mle/{project}_cnv_correction_qc_mle.bmk",
    shell:
        """
        Rscript scripts/run_crisprcleanR_qc.R
            --sgrna-summ {input.sgrna_summary} \
            --gene-summ {input.gene_summary} \
            --lib-type {params.library_type} \
            --sgrna-library {params.sgrna_library} \
            --essen-genes {params.essential_genes} \
            --non-essen-genes {params.non_essential_genes} \
            --label {params.label}
        
        """