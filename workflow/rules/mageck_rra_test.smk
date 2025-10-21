rule mageck_rra_test:
    input:
        count_table=lambda wc: get_count_table_to_test(wc.project, wc.mageck_test),
    output:
        sgrna_summary="results/mageck_rra_test/{mageck_test}/{project}_{mageck_test}.sgrna_summary.txt",
        gene_summary="results/mageck_rra_test/{mageck_test}/{project}_{mageck_test}.gene_summary.txt",
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
        extra=config["parameters"]["mageck_rra_test"]["extra"],
        out_prefix=lambda wc: f"results/mageck_rra_test/{wc.mageck_test}/{wc.project}_{wc.mageck_test}",
    log:
        "logs/mageck_rra_test/{mageck_test}/{project}_{mageck_test}_mageck_rra_test.log",
    benchmark:
        "benchmarks/mageck_rra_test/{mageck_test}/{project}_{mageck_test}_mageck_rra_test.bmk"
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


# TODO: Try to include this rule in the more general one.
rule mageck_rra_NCT_qc:
    input:
        count_table="results/filter_sgrna_counts/{project}_processed_NCT.count.txt",
    output:
        sgrna_summary="results/mageck_rra_NCT_qc/{project}_processed_NCT.sgrna_summary.txt",
        gene_summary="results/mageck_rra_NCT_qc/{project}_processed_NCT.gene_summary.txt",
        normed_count="results/mageck_rra_NCT_qc/{project}_processed_NCT.normalized.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        treat_samples=",".join(treat_samples),
        ctrl_samples=",".join(ctrl_samples),
        norm_method=config["parameters"]["mageck_rra_NCT_qc"]["norm_method"],
        gene_lfc_method=config["parameters"]["mageck_rra_test"]["gene_lfc_method"],
        gene_fdr_thres=config["parameters"]["mageck_rra_test"]["gene_fdr_thres"],
        p_adj_method=config["parameters"]["mageck_rra_test"]["p_adj_method"],
        out_prefix=lambda wc: f"results/mageck_rra_NCT_qc/{wc.project}_processed_NCT",
        extra=config["parameters"]["mageck_rra_NCT_qc"].get("extra", ""),
    log:
        "logs/mageck_rra_NCT_qc/{project}_mageck_rra_NCT_qc.log",
    benchmark:
        "benchmarks/mageck_rra_NCT_qc/{project}_mageck_rra_NCT_qc.bmk"
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
            --output-prefix {params.out_prefix} \
            {params.extra} 2> {log}
        """
