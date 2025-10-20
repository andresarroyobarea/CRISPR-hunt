rule bagel2_fc:
    input:
        processed_count = "results/filter_sgrna_counts/{project}_processed.count.txt",
    output: 
        fold_change = "results/bagel2_fc/{project}_total.foldchange",
        norm_counts = "results/bagel2_fc/{project}_total.normed_readcount",
    conda:
        config["conda_envs"]["bagel2"],
    params:
        out_prefix=lambda wc: f"results/bagel2_fc/{wc.project}",
        ctrl_samples = ",".join(ctrl_samples),
        pseudo_count = config["parameters"]["bagel2_fc"]["pseudo_count"],
        extra = config["parameters"]["bagel2_fc"]["extra"],
    resources:
        threads = get_resource("bagel2_fc", "threads"),
    log:
        "logs/bagel2/bagel2_fc/{project}_bagel2_fc.log",
    benchmark:
        "benchmarks/bagel2/bagel2_fc/{project}_bagel2_fc.bmk",
    shell:
        """
        BAGEL.py fc \
            -i {input.processed_count} \
            -o {params.out_prefix} \
            -c {params.ctrl_samples} \
            -Np {params.pseudo_count} \
            {params.extra} \
            > /dev/null 2>> {log}
        """

rule bagel2_bf:
    input:
        fold_change=lambda wc: get_lfc_file(wc.project, wc.norm_state),
        essential_genes=config["parameters"]["bagel2_bf"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_bf"]["non_essential_genes"],
    output:
        bayes_factors="results/bagel2_bf/{project}_{norm_state}.bf",
    conda:
        config["conda_envs"]["bagel2"]
    params:
        treat_samples=",".join(treat_samples),
        extra=config["parameters"]["bagel2_bf"]["extra"],
    log:
        "logs/bagel2_bf/{project}_{norm_state}_bagel2_bf.log",
    benchmark:
        "benchmarks/bagel2_bf/{project}_{norm_state}_bagel2_bf.log"
    shell:
        """
        BAGEL.py bf \
            -i {input.fold_change} \
            -o {output.bayes_factors} \
            -e {input.essential_genes} \
            -n {input.non_essential_genes} \
            -c {params.treat_samples} \
            {params.extra} 2> {log}
        """


