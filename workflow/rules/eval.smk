localrules: bgzip, tabix, reformat_truvari_results, reformat_truvari_results_svim, cat_truvari_results

def get_vcf(wildcards):
    return config["truth"][wildcards.vcf]

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"


rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


rule callset_eval_svim:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}"
    threads: 1
    log:
        "logs/truvari/truvari.svim.{data}.{aligner}.{parameters}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


# rule callset_eval_svim_gtcomp:
#     input:
#         genome = config["reference"],
#         truth_vcf = get_vcf,
#         truth_bed = config["truth"]["bed"],
#         calls = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz",
#         index = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz.tbi"
#     output:
#         temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/base-filter.vcf"),
#         temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/call-filter.vcf"),
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/fn.vcf",
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/fp.vcf",
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/giab_report.txt",
#         temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/log.txt"),
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/summary.txt",
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/tp-base.vcf",
#         "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/tp-call.vcf"
#     params:
#         out_dir="{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt"
#     threads: 1
#     log:
#         "logs/{aligner}/truvari/pooled.svim.{sample}.{parameters}.{minscore}.{vcf}.gt.log"
#     shell:
#         "rm -rf {params.out_dir} && truvari -f {input.genome}\
#                     -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
#                     --passonly --gtcomp --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule callset_eval:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        "pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}"
    threads: 1
    log:
        "logs/truvari/truvari.{caller}.{data}.{aligner}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


# rule callset_eval_gtcomp:
#     input:
#         genome = config["reference"],
#         truth_vcf = get_vcf,
#         truth_bed = config["truth"]["bed"],
#         calls = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz",
#         index = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz.tbi"
#     output:
#         temp("{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/base-filter.vcf"),
#         temp("{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/call-filter.vcf"),
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/fn.vcf",
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/fp.vcf",
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/giab_report.txt",
#         temp("{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/log.txt"),
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/summary.txt",
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/tp-base.vcf",
#         "{aligner}/{caller,Sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/tp-call.vcf"
#     params:
#         out_dir="{aligner}/{caller}_results/{sample}/{minscore}/{vcf}.gt"
#     threads: 1
#     log:
#         "logs/{aligner}/truvari/pooled.{caller}.{sample}.{minscore}.{vcf}.gt.log"
#     shell:
#         "rm -rf {params.out_dir} && truvari -f {input.genome}\
#                     -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
#                     --passonly --gtcomp --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule reformat_truvari_results:
    input:
        "pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    output:
        temp("pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/pr_rec.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/summary.txt"
    output:
        temp("pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/pr_rec.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"SVIM\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results:
    input:
        svim = expand("pipeline/SVIM_results/{{aligner}}/{data}/0.3/{minscore}/{vcf}/pr_rec.txt", data = SUBSAMPLES, minscore=[0] + list(range(1, 60, 2)), vcf=VCFS),
        sniffles = expand("pipeline/Sniffles_results/{{aligner}}/{data}/{minscore}/{vcf}/pr_rec.txt", data = SUBSAMPLES, minscore=list(range(1, 60, 2)), vcf=VCFS),
        pbsv = expand("pipeline/pbsv_results/{{aligner}}/{data}/{minscore}/{vcf}/pr_rec.txt", data = SUBSAMPLES, minscore=list(range(1, 60, 2)), vcf=VCFS)
    output:
        svim = temp("pipeline/eval/{aligner}/svim.results.txt"),
        sniffles = temp("pipeline/eval/{aligner}/sniffles.results.txt"),
        pbsv = temp("pipeline/eval/{aligner}/pbsv.results.txt"),
        all = "pipeline/eval/{aligner}/all_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

# rule plot_pr_tools:
#     input:
#         "{aligner}/eval/{sample}/{parameters}/all_results.{vcf}.txt"
#     output:
#         "{aligner}/eval/{sample}/{parameters}/tools_pr_all.{vcf}.png"
#     threads: 1
#     log:
#         "logs/{aligner}/rplot/{sample}.{parameters}.tools.pr.{vcf}.log"
#     shell:
#         "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"

# rule plot_pr_tools_multiple_coverages:
#     input:
#         results = "{aligner}/eval/{sample}/{caller}_results_multiple_coverages.{vcf}.txt",
#         coverages = "{aligner}/mosdepth/mean_coverages.txt"
#     output:
#         "{aligner}/eval/{sample}/{caller,sniffles|pbsv}_pr_multiple_coverages.{vcf}.png"
#     threads: 1
#     log:
#         "logs/{aligner}/rplot/{sample}.{caller}.pr.coverages.{vcf}.log"
#     shell:
#         "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input.results} {input.coverages} {output} > {log}"

# rule plot_pr_tools_multiple_coverages_svim:
#     input:
#         results = "{aligner}/eval/{sample}/{parameters}/svim_results_multiple_coverages.{vcf}.txt",
#         coverages = "{aligner}/mosdepth/mean_coverages.txt"
#     output:
#         "{aligner}/eval/{sample}/{parameters}/svim_pr_multiple_coverages.{vcf}.png"
#     threads: 1
#     log:
#         "logs/{aligner}/rplot/{sample}.{parameters}.svim.pr.coverages.{vcf}.log"
#     shell:
#         "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input.results} {input.coverages} {output} > {log}"
