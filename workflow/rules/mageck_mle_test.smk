rule mageck_mle_test:
    input:
        count_table = count_table, 
        design = design_matrix_filt_path,
    output: 
        sgrna_summary="results/mageck_mle_test/{project}_sgrna_summary.txt",
        gene_summary="results/mageck_mle_test/{project}_gene_summary.txt",
    conda:
        config["conda_envs"]["mageck"],
    params:
        norm_method=config["parameters"]["mageck_mle_test"]["norm_method"],
        p_adj_method=config["parameters"]["mageck_mle_test"]["p_adj_method"],
        out_prefix=config["parameters"]["mageck_mle_test"]["out_prefix"],
        samples=",".join(SAMPLES),
        extra=config["parameters"]["mageck_mle_test"]["extra"],
    resources:
        threads=get_resource("mageck_mle_test", "threads"),
    log:
        "logs/mageck_mle_test/mageck_mle_test.log",
    benchmark:
        "benchmarks/mageck_mle_test/mageck_mle_test.mbk",
    shell: 
        """
        mageck mle \
            --count-table {input.count_table} \
            --design-matrix {input.design} \
            --norm-method {params.norm_method} \
            --include-samples "{params.samples}" \
            --adjust-method {params.p_adj_method} \
            --output-prefix {params.out_prefix} \
            --threads {resources.threads} \
            {params.extra} 2> {log}
        """