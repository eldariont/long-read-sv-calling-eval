localrules: filter_svim, fix_sniffles, filter_insertions_and_deletions, filter_insertions_and_deletions_svim

# SVIM
rule run_svim:
    input:
        genome = config["reference"],
        bam="pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai="pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        "pipeline/SVIM/{aligner}/{data}/{max_distance}/variants.vcf"
    params:
        working_dir = "pipeline/SVIM/{aligner}/{data}/{max_distance}/",
        min_sv_size = config["parameters"]["min_sv_size"],
        runtime = "600",
        memory = "10000"
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
        temp("pipeline/Sniffles/{aligner}/{data}/raw_{minsupport,[0-9]+}.vcf")
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        runtime = "600",
        memory = "20000"
    threads: 3    
    conda:
        "../envs/sniffles.yaml"
    shell:
        "sniffles --mapped_reads {input.bam} --min_length {params.min_sv_size} --min_support {wildcards.minsupport} --vcf {output} --threads {threads} --genotype"

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
        vcf = "pipeline/pbsv/{aligner}/{data}/min_{minsupport,[0-9]+}.vcf",
        svsig = temp("pipeline/pbsv/{aligner}/{data}/min_{minsupport}.svsig.gz")
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        runtime = "600",
        memory = "20000"
    threads: 2
    conda:
        "../envs/pbsv.yaml"
    shell:
        """
        pbsv discover {input.bam} {output.svsig} && \
        pbsv call -t INS,DEL \
        -j {threads} \
        --min-sv-length {params.min_sv_size} \
        --max-ins-length 100K \
        --call-min-reads-one-sample {wildcards.minsupport} \
        --call-min-reads-all-samples {wildcards.minsupport} \
        --call-min-reads-per-strand-all-samples 0 \
        --call-min-bnd-reads-all-samples 0 \
        --call-min-read-perc-one-sample 0 {input.genome} {output.svsig} {output.vcf}
        """

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
