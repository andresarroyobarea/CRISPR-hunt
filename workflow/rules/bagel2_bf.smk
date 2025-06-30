rule bagel2_bf:
    input:
        fold_change="results/bagel2_fc/{project}.foldchange",
    output: 
        bayes_factors="results/bagel2_bf/{project}.bf",
    conda:
        config["conda_envs"]["bagel2"],
    params:
        essential_genes=config["parameters"]["bagel2_bf"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_bf"]["non_essential_genes"],
        treat_samples=",".join(treat_samples),
        extra=config["parameters"]["bagel2_bf"]["extra"],
    log:
        "logs/bagel2_bf/{project}_bagel2_bf.log", 
    benchmark:
        "benchmarks/bagel2_bf/{project}_bagel2_bf.log",
    shell:
        """
        BAGEL.py bf \
            -i {input.fold_change} \
            -o {output.bayes_factors} \
            -e {params.essential_genes} \
            -n {params.non_essential_genes} \
            -c {params.treat_samples} \
            {params.extra} 2> {log}
        """