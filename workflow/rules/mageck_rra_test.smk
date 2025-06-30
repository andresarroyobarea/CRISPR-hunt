rule mageck_rra_test:
    input:
        count_table,
    output:
        sgrna_summary="results/mageck_rra_test/{project}_sgrna_summary.txt",
        gene_summary="results/mageck_rra_test/{project}_gene_summary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        treat_samples=",".join(treat_samples),
        ctrl_samples=",".join(ctrl_samples),
        norm_method=config["parameters"]["mageck_rra_test"]["norm_method"],
        gene_lfc_method=config["parameters"]["mageck_rra_test"]["gene_lfc_method"],
        p_adj_method=config["parameters"]["mageck_rra_test"]["p_adj_method"],
        gene_fdr_thres=config["parameters"]["mageck_rra_test"]["gene_fdr_thres"],
        remove_zeros=config["parameters"]["mageck_rra_test"]["remove_zero"],
        out_prefix=config["parameters"]["mageck_rra_test"]["out_prefix"],
        extra=config["parameters"]["mageck_rra_test"]["extra"],
    log:
        "logs/mageck_rra_test/mageck_rra_test.log",
    benchmark:
        "benchmarks/mageck_rra_test/mageck_rra_test.bmk"
    shell:
        """
        mageck test \
            --count-table {input.count_table} \
            --treatment-id {params.treat_samples} \
            --control-id {params.ctrl_samples} \
            --norm-method {params.norm_method} \
            --gene-lfc-method {params.gene_lfc_method} \
            --adjust-method {params.p_adj_method} \
            --gene-test-fdr-threshold {params.gene_fdr_thres} \
            --remove-zero {params.remove_zeros} \
            --output-prefix {params.out_prefix} \
            {params.extra} 2> {log}
        """
