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
        --output-prefix \
        {extra} 2> {log}
        """
