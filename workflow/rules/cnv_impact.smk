rule cnv_impact:
    input:
        mageck_corr="results/mageck_rra_test/{project}_norm_cnv_correction_gene_summary.txt",
        mageck_uncorr="results/mageck_rra_test/{project}_{mageck_norms}_gene_summary.txt",
    output:
        cnv_impact_stats="results/cnv_impact/{mageck_norms}/{project}_{mageck_norms}_cnv_impact_stats.csv",
        cnv_impact_plots="results/cnv_impact/{mageck_norms}/{project}_{mageck_norms}_cnv_impact_plots.pdf",
    conda:
        config["conda_envs"]["crisprcleanR"]
    params:
        fdr_thres=config["parameters"]["cnv_impact"]["fdr_cnv"],
        label=config["project"],
    log:
        "logs/cnv_impact/{project}_{mageck_norms}_cnv_correction_impact.log",
    benchmark:
        "benchmarks/cnv_impact/{project}_{mageck_norms}_cnv_correction_impact.bmk"
    shell:
        """
        Rscript scripts/run_cnv_impact.R \
            --mageck-gene-summ-corr {input.mageck_corr} \
            --mageck-gene-summ-uncorr {input.mageck_uncorr} \
            --fdr-threshold {params.fdr_cnv} \
            --label {params.label} \
            --out-cnv-impact-stats {output.cnv_impact_stats} \
            --out-cnv-impact-plots {output.cnv_impact_plots} 2> {log}
        """
