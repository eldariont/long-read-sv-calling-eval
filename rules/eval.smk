def get_vcf(wildcards):
    return config["truth"][wildcards.vcf]

rule sort_calls:
    input:
        "{aligner}/{caller}_calls/{sample}.min_{minscore}.vcf"
    output:
        temp("{aligner}/{caller,sniffles|pbsv}_calls/{sample}.min_{minscore,[0-9]+}.sorted.vcf")
    threads: 1
    log:
        "logs/{aligner}/bcftools_sort/sorting_{caller}_{sample}_{minscore}.log"
    shell:
        "bcftools sort {input} > {output} 2> {log}"

rule sort_calls_svim:
    input:
        "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.vcf"
    output:
        temp("{aligner}/svim_calls/{sample}/{parameters}/min_{minscore,[0-9]+}.truvari.sorted.vcf")
    threads: 1
    log:
        "logs/{aligner}/bcftools_sort/sorting_svim_{sample}_{parameters}_{minscore}.truvari.log"
    shell:
        "bcftools sort {input} > {output} 2> {log}"

rule bgzip:
    input:
        "{name}.sorted.vcf"
    output:
        "{name}.sorted.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"


rule tabix:
    input:
        "{name}.sorted.vcf.gz"
    output:
        "{name}.sorted.vcf.gz.tbi"
    shell:
        "tabix {input}"


rule callset_eval_svim:
    input:
        genome = config["genome"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz",
        index = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz.tbi"
    output:
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/base-filter.vcf"),
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/call-filter.vcf"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/fn.vcf",
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/fp.vcf",
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/giab_report.txt",
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/log.txt"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/summary.txt",
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/tp-base.vcf"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/tp-call.vcf"
    params:
        out_dir="{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.svim.{sample}.{parameters}.{minscore}.{vcf}.log"
    shell:
        "rm -rf {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule callset_eval_svim_gtcomp:
    input:
        genome = config["genome"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz",
        index = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz.tbi"
    output:
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/base-filter.vcf"),
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/call-filter.vcf"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/fn.vcf",
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/fp.vcf",
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/giab_report.txt",
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/log.txt"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/summary.txt",
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/tp-base.vcf"),
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt/tp-call.vcf"
    params:
        out_dir="{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}.gt"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.svim.{sample}.{parameters}.{minscore}.{vcf}.gt.log"
    shell:
        "rm -rf {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --gtcomp --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule callset_eval:
    input:
        genome = config["genome"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/{caller}_calls/{sample}.min_{minscore}.sorted.vcf.gz",
        index = "{aligner}/{caller}_calls/{sample}.min_{minscore}.sorted.vcf.gz.tbi"
    output:
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/base-filter.vcf"),
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/call-filter.vcf"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/fn.vcf",
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/fp.vcf",
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/giab_report.txt",
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/log.txt"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/summary.txt",
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/tp-base.vcf"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}/tp-call.vcf"
    params:
        out_dir="{aligner}/{caller}_results/{sample}/{minscore}/{vcf}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.{caller}.{sample}.{minscore}.{vcf}.log"
    shell:
        "rm -rf {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule callset_eval_gtcomp:
    input:
        genome = config["genome"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/{caller}_calls/{sample}.min_{minscore}.sorted.vcf.gz",
        index = "{aligner}/{caller}_calls/{sample}.min_{minscore}.sorted.vcf.gz.tbi"
    output:
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/base-filter.vcf"),
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/call-filter.vcf"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/fn.vcf",
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/fp.vcf",
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/giab_report.txt",
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/log.txt"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/summary.txt",
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/tp-base.vcf"),
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/{vcf}.gt/tp-call.vcf"
    params:
        out_dir="{aligner}/{caller}_results/{sample}/{minscore}/{vcf}.gt"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.{caller}.{sample}.{minscore}.{vcf}.gt.log"
    shell:
        "rm -rf {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --gtcomp --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"


rule reformat_truvari_results:
    input:
        "{aligner}/{caller}_results/{sample}/{minscore}/{vcf}/summary.txt"
    output:
        temp("{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/pr_rec.{vcf}.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.sample}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/{vcf}/summary.txt"
    output:
        temp("{aligner}/svim_results/{sample}/{parameters}/{minscore}/pr_rec.{vcf}.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"svim\", \"{wildcards.sample}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results_all:
    input:
        expand("{{aligner}}/svim_results/{{sample}}/{{parameters}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("{{aligner}}/sniffles_results/{{sample}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("{{aligner}}/pbsv_results/{{sample}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"]))
    output:
        "{aligner}/eval/{sample}/{parameters}/all_results.{vcf}.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_svim_multiple_coverages:
    input:
        expand("{{aligner}}/svim_results/{{sample}}.subsampled.{fraction}/{{parameters}}/{minscore}/pr_rec.{{vcf}}.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("{{aligner}}/svim_results/{{sample}}/{{parameters}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"]))
    output:
        "{aligner}/eval/{sample}/{parameters}/svim_results_multiple_coverages.{vcf}.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_sniffles_multiple_coverages:
    input:
        expand("{{aligner}}/sniffles_results/{{sample}}.subsampled.{fraction}/{minscore}/pr_rec.{{vcf}}.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("{{aligner}}/sniffles_results/{{sample}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"]))
    output:
        "{aligner}/eval/{sample}/sniffles_results_multiple_coverages.{vcf}.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_pbsv_multiple_coverages:
    input:
        expand("{{aligner}}/pbsv_results/{{sample}}.subsampled.{fraction}/{minscore}/pr_rec.{{vcf}}.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"])),
        expand("{{aligner}}/pbsv_results/{{sample}}/{minscore}/pr_rec.{{vcf}}.txt", minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"]))
    output:
        "{aligner}/eval/{sample}/pbsv_results_multiple_coverages.{vcf}.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule plot_pr_tools:
    input:
        "{aligner}/eval/{sample}/{parameters}/all_results.{vcf}.txt"
    output:
        "{aligner}/eval/{sample}/{parameters}/tools_pr_all.{vcf}.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{parameters}.tools.pr.{vcf}.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"

rule plot_pr_tools_multiple_coverages:
    input:
        results = "{aligner}/eval/{sample}/{caller}_results_multiple_coverages.{vcf}.txt",
        coverages = "{aligner}/mosdepth/mean_coverages.txt"
    output:
        "{aligner}/eval/{sample}/{caller,sniffles|pbsv}_pr_multiple_coverages.{vcf}.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{caller}.pr.coverages.{vcf}.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input.results} {input.coverages} {output} > {log}"

rule plot_pr_tools_multiple_coverages_svim:
    input:
        results = "{aligner}/eval/{sample}/{parameters}/svim_results_multiple_coverages.{vcf}.txt",
        coverages = "{aligner}/mosdepth/mean_coverages.txt"
    output:
        "{aligner}/eval/{sample}/{parameters}/svim_pr_multiple_coverages.{vcf}.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{parameters}.svim.pr.coverages.{vcf}.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input.results} {input.coverages} {output} > {log}"
