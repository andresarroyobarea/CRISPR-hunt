rule bagel2_fc:
    input:
        count_table_raw="results/mageck_count_raw/raw.counts.txt",
    output:
        fold_change="results/bagel2_fc/{project}.foldchange",
        norm_counts="results/bagel2_fc/{project}.normed_readcount",
    conda:
        config["conda_envs"]["bagel2"]
    params:
        project=config["project"],
        ctrl_samples=",".join(ctrl_samples),
        extra=config["parameters"]["bagel2_fc"]["extra"],
    log:
        "logs/bagel2_fc/{project}_bagel2_fc.log",
    benchmark:
        "benchmarks/bagel2/{project}_bagel2_fc.bmk"
    shell:
        """
        BAGEL.py fc \
            -i {input.count_table_raw} \
            -o {params.project} \
            -c {params.ctrl_samples} \
            {params.extra} 2> {log}
        """
