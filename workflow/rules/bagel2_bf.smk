rule bagel2_bf:
    input:
        fold_change="results/bagel2_fc/{project}.foldchange",
        essential_genes=config["parameters"]["bagel2_bf"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_bf"]["non_essential_genes"],
    output:
        bayes_factors="results/bagel2_bf/{project}.bf",
    conda:
        config["conda_envs"]["bagel2"]
    params:
        treat_samples=",".join(treat_samples),
        extra=config["parameters"]["bagel2_bf"]["extra"],
    log:
        "logs/bagel2_bf/{project}_bagel2_bf.log",
    benchmark:
        "benchmarks/bagel2_bf/{project}_bagel2_bf.log"
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
