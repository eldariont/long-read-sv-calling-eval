localrules: filter_svim, fix_sniffles, filter_insertions_and_deletions, filter_insertions_and_deletions_svim

# SVIM
rule run_svim:
    input:
        genome = config["reference"],
        bam="pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai="pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        "pipeline/SVIM/{aligner}/{data}/{max_distance}/variants.vcf"
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "pipeline/SVIM/{aligner}/{data}/{max_distance}/",
        min_sv_size = config["parameters"]["min_sv_size"]
    threads: 1
    shell:
        "/home/heller_d/bin/anaconda3/bin/svim alignment --sample {wildcards.data} --cluster_max_distance {wildcards.max_distance} \
         --min_sv_size {params.min_sv_size} --segment_gap_tolerance 20 --segment_overlap_tolerance 20 \
         --interspersed_duplications_as_insertions --tandem_duplications_as_insertions --read_names {params.working_dir} {input.bam} {input.genome}"

rule filter_svim:
    input:
        "pipeline/SVIM/{aligner}/{data}/{max_distance}/variants.vcf"
    output:
        temp("pipeline/SVIM/{aligner}/{data}/{max_distance}/min_{minscore,[0-9]+}.vcf")
    threads: 1
    shell:
        "grep -v \"hom_ref\" {input} | \
         awk 'OFS=\"\\t\" {{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>={wildcards.minscore}) {{ print $1, $2, $3, $4, $5, $6, \"PASS\", $8, $9, $10 }} }} }}' > {output}"

# SNIFFLES
rule run_sniffles:
    input:
        bam = "pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai = "pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        expand("pipeline/Sniffles/{{aligner}}/{{data}}/raw_{minsupport}.vcf", 
                minsupport=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])))
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        tmpdir = "300",
        sniffles_from = config["minimums"]["sniffles_from"],
        sniffles_to = config["minimums"]["sniffles_to"],
        sniffles_step = config["minimums"]["sniffles_step"],
        outdir = "pipeline/Sniffles/{aligner}/{data}/"
    threads: 30
    conda:
        "../envs/sniffles.yaml"
    shell:
        "bash workflow/scripts/run_sniffles.sh {input.bam} {input.bai} {params.sniffles_from} {params.sniffles_to} {params.sniffles_step} {params.min_sv_size} {threads} {params.outdir}"

#see https://github.com/spiralgenetics/truvari/issues/43
rule fix_sniffles:
    input:
        "pipeline/Sniffles/{aligner}/{data}/raw_{support}.vcf"
    output:
        "pipeline/Sniffles/{aligner}/{data}/min_{support,[0-9]+}.vcf"
    shell:
        "sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}"

#PBSV
rule run_pbsv:
    input:
        bam = "pipeline/alignment_pooled/{data}.pbmm2.bam",
        bai = "pipeline/alignment_pooled/{data}.pbmm2.bam.bai",
        genome = config["reference"],
    output:
        expand("pipeline/pbsv/{{aligner}}/{{data}}/min_{minsupport}.vcf",
                minsupport=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])))
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        tmpdir = "500",
        pbsv_from = config["minimums"]["pbsv_from"],
        pbsv_to = config["minimums"]["pbsv_to"],
        pbsv_step = config["minimums"]["pbsv_step"],
        outdir = "pipeline/pbsv/{aligner}/{data}/"
    threads: 15
    conda:
        "../envs/pbsv.yaml"
    shell:
        "bash workflow/scripts/run_pbsv.sh {input.bam} {input.bai} {input.genome} {params.pbsv_from} {params.pbsv_to} {params.pbsv_step} {params.min_sv_size} {threads} {params.outdir}"

#Split to SV classes
rule filter_insertions_and_deletions:
    input:
        "pipeline/{caller}/{aligner}/{data}/min_{minscore}.vcf"
    output:
        "pipeline/{caller,Sniffles|pbsv}/{aligner}/{data}/min_{minscore,[0-9]+}.indel.vcf"
    threads: 1
    shell:
        "bcftools sort {input} | bcftools view -i 'SVTYPE=\"DEL\" | SVTYPE=\"INS\"' -Ov > {output}"

rule filter_insertions_and_deletions_svim:
    input:
        "pipeline/SVIM/{aligner}/{data}/{max_distance}/min_{minscore}.vcf"
    output:
        "pipeline/SVIM/{aligner}/{data}/{max_distance}/min_{minscore,[0-9]+}.indel.vcf"
    threads: 1
    shell:
        "bcftools sort {input} | bcftools view -i 'SVTYPE=\"DEL\" | SVTYPE=\"INS\"' -Ov > {output}"
