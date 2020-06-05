localrules: plot_pr_all_results, plot_pr_tools, plot_pr_coverages, plot_pr_svim_parameters, SV_length_plot, SV_length_plot_svim

rule plot_pr_all_results:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.all.png"
    threads: 1
    log:
        "pipeline/logs/rplot.all.{aligner}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_all_results.R {input} {output} > {log}"

rule plot_pr_tools:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.tools.{vcf}.png"
    threads: 1
    log:
        "pipeline/logs/rplot.tools.{aligner}.{vcf}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools.R {input} {wildcards.vcf} {output} > {log}"

rule plot_pr_coverages:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        png = "pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.png",
        tsv = "pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.tsv"
    threads: 1
    log:
        "pipeline/logs/rplot.coverages.{aligner}.{vcf}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_coverages.R {input} {wildcards.vcf} {output.png} {output.tsv} > {log}"

rule plot_pr_svim_parameters:
    input:
        "pipeline/eval/{aligner}/svim_parameter_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.svim.parameters.png"
    threads: 1
    log:
        "pipeline/logs/rplot.all.{aligner}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_svim_parameters.R {input} {output} > {log}"

# rule plot_pr_coverages_bar:
#     input:
#         "pipeline/eval/{aligner}/all_results.txt"
#     output:
#         "pipeline/eval/{aligner}/results.{aligner}.coverages.bar.png"
#     threads: 1
#     log:
#         "pipeline/logs/rplot.coveragesbar.{aligner}.log"
#     shell:
#         "Rscript --vanilla workflow/scripts/plot_coverages_bar.R {input} {output} > {log}"

rule SV_length_plot:
    input:
        "pipeline/{caller}/{aligner}/{data}/min_{minimum}.indel.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_{caller,Sniffles|pbsv}_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-nucleotides_affected_{caller,sniffles|pbsv}_{minimum}.txt",
    log:
        "logs/svplot/svlength_{caller}_{aligner}_{data}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"

rule SV_length_plot_svim:
    input:
        "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minimum}.indel.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_SVIM_{parameters}_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-nucleotides_affected_SVIM_{parameters}_{minimum}.txt",
    log:
        "logs/svplot/svlength_SVIM_{aligner}_{data}_{parameters}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"