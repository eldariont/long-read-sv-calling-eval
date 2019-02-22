def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule minimap2_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        "minimap2/alignment/{sample}.bam"
    params:
        preset = config["parameters"]["minimap_preset"]
    threads:
        8
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 --MD -ax {params.preset} -t {threads} {input.genome} {input.fq} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule pbmm_index:
    input:
        genome = config["genome"]
    output:
        index = config["genome"] + ".mmi"
    params:
        preset = config["parameters"]["pbmm_preset"]
    threads: 1
    shell:
        "pbmm2 index {input.genome} {output.index} --preset {params.preset}"

rule pbmm_align:
    input:
        fq = get_samples,
        index = config["genome"] + ".mmi"
    output:
        "minimap2_pbsv/alignment/{sample}.bam"
    threads:
        8
    params:
        sample = "{sample}",
        preset = config["parameters"]["pbmm_preset"]
    log:
        "logs/minimap2_pbsv/{sample}.log"
    shell:
        """
        pbmm2 align {input.index} {input.fq} {output} --preset {params.preset} \
                --sort --rg '@RG\tID:rg1a\tSM:{params.sample}' --sample HG2 2> {log}
        """

rule ngmlr_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        protected("ngmlr/alignment/{sample}.bam")
    params:
        preset = config["parameters"]["ngmlr_preset"]
    threads:
        36
    log:
        "logs/ngmlr/{sample}.log"
    shell:
        "zcat {input.fq} | \
         ngmlr --presets {params.preset} -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule minimap2_last_like_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        "minimap2_last_like/alignment/{sample}.bam"
    params:
        preset = config["parameters"]["minimap_preset"]
    threads:
        8
    log:
        "logs/minimap2_last_like/{sample}.log"
    shell:
        "minimap2 --MD -ax {params.preset} --no-long-join -r50 -t {threads} {input.genome} {input.fq} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule samtools_index:
    input:
        "{aligner}/alignment/{sample}.bam"
    output:
        "{aligner}/alignment/{sample}.bam.bai"
    threads: 4
    log:
        "logs/{aligner}/samtools_index/{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"

rule alignment_stats:
    input:
        bam = expand("{{aligner}}/alignment/{sample}.bam", sample=config["samples"]),
        bai = expand("{{aligner}}/alignment/{sample}.bam.bai", sample=config["samples"])
    output:
        "{aligner}/alignment_stats/alignment_stats.txt"
    log:
        "logs/{aligner}/alignment_stats/alignment_stats.log"
    shell:
        os.path.join(workflow.basedir, "scripts/alignment_stats.py") + \
            " -o {output} {input.bam} 2> {log}"

rule pool_samples:
    input:
        expand("{{aligner}}/alignment/{sample}.bam", sample=config["samples"])
    output:
        "{aligner}/alignment_pooled/pooled.bam"
    log:
        "logs/{aligner}/pool_samples.log"
    shell:
        "samtools merge -r {output} {input}"

rule subsample_alignments:
    input:
        "{aligner}/alignment_pooled/pooled.bam"
    output:
        "{aligner}/alignment_pooled/pooled.subsampled.{fraction}.bam"
    threads: 4
    params:
        additional_threads = 3
    log:
        "logs/{aligner}/samtools_view/subsample.{fraction}.log"
    shell:
        "samtools view -s 10.{wildcards.fraction} -@ {params.additional_threads} -b {input} -o {output}"
