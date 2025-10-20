rule bagel2_fc:
    input:
        processed_count = "results/filter_sgrna_counts/{project}_processed.count.txt",
    output: 
        fold_change = "results/bagel2_fc/{project}.foldchange",
        norm_counts = "results/bagel2_fc/{project}.normed_readcount",
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
