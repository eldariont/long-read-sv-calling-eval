rule sort_calls:
    input:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.vcf"
    output:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf"
    threads: 1
    log:
        "logs/{aligner}/bcftools_sort/sorting_{caller}_{minscore}.log"
    shell:
        "bcftools sort {input} > {output} 2> {log}"

rule bgzip:
    input:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf"
    output:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"


rule tabix:
    input:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf.gz"
    output:
        "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf.gz.tbi"
    shell:
        "tabix {input}"


rule callset_eval:
    input:
        genome = config["genome"],
        truth = config["truth"],
        calls = "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf.gz",
        index = "{aligner}/{caller}_calls/pooled.min_{minscore}.sorted.vcf.gz.tbi"
    output:
        "{aligner}/{caller}_results/{minscore}/summary.txt"
    params:
        out_dir="{aligner}/{caller}_results/{minscore}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.{caller}.{minscore}.log"
    shell:
        "rm -r {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth} -c {input.calls} -o {params.out_dir}\
                    --passonly -r 1000 -p 0.00 2> {log}"


rule reformat_truvari_results:
    input:
        "{aligner}/{caller}_results/{minscore}/summary.txt"
    output:
        "{aligner}/{caller}_results/{minscore}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat truvari_results:
    input:
        expand("{{aligner}}/svim_results/{minscore}/pr_rec.txt", minscore=range(1, 100, 5)),
        expand("{{aligner}}/sniffles_results/{minscore}/pr_rec.txt", minscore=range(1, 42, 5))
    output:
        "{aligner}/eval/all_results.txt"
    threads: 1
    shell:
        "cat {input} > {output}"


rule plot_pr_tools:
    input:
        "{aligner}/eval/all_results.txt"
    output:
        "{aligner}/eval/tools_pr.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/pooled.tools.pr.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"
