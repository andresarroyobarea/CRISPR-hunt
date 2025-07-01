rule bagel2_pr:
    input:
        bayes_factors="results/bagel2_bf/{project}.bf",
        essential_genes=config["parameters"]["bagel2_bf"]["essential_genes"],
        non_essential_genes=config["parameters"]["bagel2_bf"]["non_essential_genes"],
    output:
        prec_recall="results/bagel2_pr/{project}-pr",
    conda:
        config["conda_envs"]["bagel2"]
    params:
        extra=config["parameters"]["bagel2_pr"]["extra"],
    log:
        "logs/bagel2_pr/{project}_bagel2_pr.log",
    benchmark:
        "benchmarks/bagel2_pr/{project}_bagel2_pr.bmk"
    shell:
        """
        BAGEL.py pr \
            -i {input.bayes_factors} \
            -o {output.prec_recall} \
            -e {input.essential_genes} \
            -n {input.non_essential_genes} \
            {params.extra} 2> {log}
        """
