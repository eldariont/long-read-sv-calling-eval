configfile: "config.yaml"

rule mosdepth_get:
    input:
        bam = "{aligner}/alignment_pooled/{sample}.bam",
        bai = "{aligner}/alignment_pooled/{sample}.bam.bai"
    threads: 4
    output:
        protected("{aligner}/mosdepth/{sample}.mosdepth.global.dist.txt"),
        protected("{aligner}/mosdepth/{sample}.regions.bed.gz"),
    params:
        windowsize = 500,
        prefix = "{sample}",
        aligner = "{aligner}"
    log:
        "logs/{aligner}/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  {params.aligner}/mosdepth/{params.prefix} {input.bam} 2> {log}"


rule mosdepth_combine:
    input:
        expand("{{aligner}}/mosdepth/{subset}.regions.bed.gz", subset=["pooled.subsampled.{0}".format(cov) for cov in range(10, 91, 10)] + ["pooled"])
    output:
        "{aligner}/mosdepth/regions.combined.gz"
    log:
        "logs/{aligner}/mosdepth/mosdepth_combine.log"
    shell:
        os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
            " {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:
        expand("{{aligner}}/mosdepth/{subset}.mosdepth.global.dist.txt", subset=["pooled.subsampled.{0}".format(cov) for cov in range(10, 91, 10)] + ["pooled"])
    output:
        "{aligner}/mosdepth_global_plot/global.html"
    log:
        "logs/{aligner}/mosdepth/mosdepth_global_plot.log"
    shell:
        os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
            " {input} -o {output} 2> {log}"


rule compute_mean_converage:
    input:
        expand("{{aligner}}/mosdepth/{subset}.mosdepth.global.dist.txt", subset=["pooled.subsampled.{0}".format(cov) for cov in range(10, 91, 10)] + ["pooled"])
    output:
        "{aligner}/mosdepth/mean_coverages.txt"
    log:
        "logs/{aligner}/mosdepth/mosdepth_mean_coverage.log"
    run:
        shell("rm -f {output}")
        for f in input:
            short_name = f.split("/")[2][:-25]
            shell("echo -n -e \"{short_name}\t\" >> {output}")
            shell("python scripts/mosdepth_mean-depth.py {f} >> {output}")