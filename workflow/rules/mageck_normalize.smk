rule mageck_normalize:
    input:
        counts_filt="results/filter_sgrna_counts/{project}_processed.count.txt",
    output:
        count_norm="results/mageck_normalize/{mageck_norms}/{project}_{mageck_norms}.count_normalized.txt",
        count_summary="results/mageck_normalize/{mageck_norms}/{project}_{mageck_norms}.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        norm_method=lambda wc: config["mageck_normalize"][wc.mageck_norms][
            "norm_method"
        ],
        extra=lambda wc: config["mageck_normalize"][wc.mageck_norms]["extra"],
        out_prefix=lambda wc: f"results/mageck_normalize/{wc.mageck_norms}/{wc.project}_{wc.mageck_norms}",
    log:
        "logs/mageck_normalize/{mageck_norms}_mageck_normalize.log",
    benchmark:
        "benchmarks/mageck_normalize/{mageck_norms}_mageck_normalize.bmk"
    shell:
        """
        mageck count \
            --count-table {output.counts_filt} \
            --norm-method {params.norm_method} \
            --output-prefix {params.out_prefix} \
            {params.extra} 2> {log}
        """
