
rule time_svim:
    input:
        genome = config["reference"],
        bam="pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai="pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        "pipeline/runtimes/SVIM/{aligner}/{data}/{pmd,[0-9]+}_{dn,[0-9]+}_{cmd,[0-9\.]+}_run{run}/SVIM.times"
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "pipeline/runtimes/SVIM/{aligner}/{data}/{pmd}_{dn}_{cmd}_run{run}/",
        min_sv_size = config["parameters"]["min_sv_size"]
    threads: 1
    shell:
        "(/usr/bin/time -v sh -c '/home/heller_d/bin/anaconda3/bin/svim alignment --sample {wildcards.data} \
         --partition_max_distance {wildcards.pmd} \
         --distance_normalizer {wildcards.dn} \
         --cluster_max_distance {wildcards.cmd} \
         --min_sv_size {params.min_sv_size} \
         --segment_gap_tolerance 20 \
         --segment_overlap_tolerance 20 \
         --interspersed_duplications_as_insertions \
         --tandem_duplications_as_insertions \
         --read_names \
         --verbose \
         {params.working_dir} {input.bam} {input.genome}') 2>> {output}"


rule time_pbsv:
    input:
        bam = "pipeline/alignment_pooled/{data}.pbmm2.bam",
        bai = "pipeline/alignment_pooled/{data}.pbmm2.bam.bai",
        genome = config["reference"],
    output:
        svsig = temp("pipeline/runtimes/pbsv/{aligner}/{data}_run{run}/svsig_all_types.svsig.gz"),
        vcf = "pipeline/runtimes/pbsv/{aligner}/{data}_run{run}/min_5.vcf",
        time = "pipeline/runtimes/pbsv/{aligner}/{data}_run{run}/pbsv.times"
    params:
        min_sv_size = config["parameters"]["min_sv_size"]
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    threads: 1
    conda:
        "../envs/pbsv.yaml"
    shell:
        """
        (/usr/bin/time -v sh -c 'pbsv discover {input.bam} {output.svsig} && pbsv call -t INS,DEL -j 1 \
        --min-sv-length {params.min_sv_size} \
        --max-ins-length 100K \
        --call-min-reads-one-sample 5 \
        --call-min-reads-all-samples 5 \
        --call-min-reads-per-strand-all-samples 0 \
        --call-min-bnd-reads-all-samples 0 \
        --call-min-read-perc-one-sample 0 {input.genome} {output.svsig} {output.vcf}') 2>> {output.time}
        """

rule time_sniffles:
    input:
        bam = "pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai = "pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        vcf = "pipeline/runtimes/Sniffles/{aligner}/{data}_run{run}/raw_5.vcf",
        time = "pipeline/runtimes/Sniffles/{aligner}/{data}_run{run}/Sniffles.times"
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
    threads: 1
    conda:
        "../envs/sniffles.yaml"
    shell:
        "(/usr/bin/time -v sh -c 'sniffles --mapped_reads {input.bam} --min_length {params.min_sv_size} --min_support 5 --vcf {output.vcf} --threads 1 --genotype') 2>> {output.time}"


