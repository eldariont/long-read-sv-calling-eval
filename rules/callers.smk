def get_chromosomes(genome):
    '''
    Gets the chromosome identifiers from the fasta genome
    '''
    fai = genome + ".fai"
    if not os.path.isfile(fai):
        sys.exit("Fasta index {} not found".format(fai))
    fa_chr = [i.split('\t')[0] for i in open(fai)]
    return fa_chr


CHROMOSOMES = get_chromosomes(config["genome"])

rule svim_call:
    input:
        "{aligner}/alignment_pooled/{sample}.bam"
    output:
        "{aligner}/svim_calls/{sample}/final_results.vcf"
    threads: 1
    log:
        "logs/{aligner}/svim_call/{sample}.log"
    shell:
        "svim alignment --sample {wildcards.sample} \
         {wildcards.aligner}/svim_calls/{wildcards.sample}/ {input} 2> {log}"

rule filter_svim:
    input:
        "{aligner}/svim_calls/{sample}/final_results.vcf"
    output:
        "{aligner}/svim_calls/{sample}.min_{minscore,[0-9]+}.vcf"
    threads: 1
    log:
        "logs/{aligner}/svim_call/{sample}.filter.{minscore}.log"
    shell:
        "cat {input} | \
         awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>={wildcards.minscore}) {{ print $0 }} }} }}' > {output}"

rule reformat_svim_calls_for_truvari:
    input:
        "{aligner}/svim_calls/{sample}.min_{minscore,[0-9]+}.vcf"
    output:
        "{aligner}/svim_calls/{sample}.min_{minscore,[0-9]+}.truvari.vcf"
    threads: 1
    shell:
        "cat {input} | sed 's/INS:NOVEL/INS/g' | sed 's/DUP:INT/INS/g' | sed 's/DUP:TANDEM/INS/g' | \
         awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($5==\"<DEL>\" || $5==\"<INS>\") {{ print $0 }} }} }}' > {output}"

rule sniffles_call:
    input:
        "{aligner}/alignment_pooled/{sample}.bam"
    output:
        "{aligner}/sniffles_calls/{sample}.min_{minsupport,[0-9]+}.vcf"
    threads: 1
    log:
        "logs/{aligner}/sniffles_call/{sample}.{minsupport}.log"
    shell:
        "sniffles --mapped_reads {input} --min_support {wildcards.minsupport} --vcf {output} --threads {threads} > {log}"

rule samtools_split:
    input:
        bam = "{aligner}/alignment_pooled/{sample}.bam",
        bai = "{aligner}/alignment_pooled/{sample}.bam.bai",
    output:
        temp("{aligner}/alignment_pooled/{sample}-{chromosome}.bam")
    params:
        chrom = "{chromosome}"
    log:
        "logs/{aligner}/samtools_split/{sample}-{chromosome}.log"
    shell:
        "samtools view {input.bam} {params.chrom} -o {output} 2> {log}"

rule nanosv_call:
    '''

    call variants using NanoSV on separate chromosomes
    the shell command will first check if there are reads in this chromosome
    and if not, will just touch the output and leave it empty
    without raising an error
    '''
    input:
        bam = "{aligner}/alignment_pooled/{sample}-{chromosome}.bam",
        bai = "{aligner}/alignment_pooled/{sample}-{chromosome}.bam.bai"
    output:
        temp("{aligner}/split_nanosv_calls/{sample}-{chromosome}.vcf")
    params:
        samtools = "samtools"
    threads:
        2
    log:
        "logs/{aligner}/nanosv/{sample}-{chromosome}.log"
    shell:
        """
        reads=$(samtools idxstats {input.bam} | \
          awk 'BEGIN {{FS = "\\t"}} ; {{sum+=$3}} END {{print sum}}')
        if [ "$reads" -eq "0" ]; then
            echo "##fileformat=VCFv4.1" > {output} && \
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> {output}
            echo "NanoSV: No reads in {input.bam}" >> exceptions.txt 2> {log}
        else
            NanoSV --threads {threads} \
                    --sambamba {params.samtools} {input.bam} \
                    -o {output} 2> {log}
        fi
        """

rule nanosv_reheader:
    input:
        "{aligner}/split_nanosv_calls/{sample}-{chromosome}.vcf",
    output:
        vcf = temp("{aligner}/split_nanosv_calls_renamed/{sample}-{chromosome}.vcf"),
        sample = temp("{aligner}/split_nanosv_calls_renamed/sample_{sample}-{chromosome}.txt")
    params:
        sample = "{sample}"
    log:
        "logs/{aligner}/bcftools_reheader/{sample}-{chromosome}.log"
    shell:
        """
        echo {params.sample} > {output.sample} &&
        bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
        """

rule nanosv_cat:
    input:
        expand("{{aligner}}/split_nanosv_calls_renamed/{{sample}}-{chromosome}.vcf",
               chromosome=CHROMOSOMES)
    output:
        "{aligner}/nanosv_calls/{sample}.vcf"
    log:
        "logs/{aligner}/bcftools-concat/{sample}.log"
    shell:
        "bcftools concat {input} | bcftools sort - -o {output} 2> {log}"

# rule npinv:
    # input:
        # bam = "{aligner}/alignment_pooled/{sample}.bam",
        # bai = "{aligner}/alignment_pooled/{sample}.bam.bai",
    # output:
        # "{aligner}/npinv/{sample}.vcf"
    # log:
        # "logs/{aligner}/npinv/{sample}.log"
    # shell:
        # "npinv --input {input} --output {output}"

# rule pbsv:
    # input:
        # bam = "minimap2_pbsv/alignment_pooled/{sample}.bam",
        # bai = "minimap2_pbsv/alignment_pooled/{sample}.bam.bai",
        # genome = config["genome"],
    # output:
        # vcf = "minimap2_pbsv/pbsv/{sample}.vcf",
        # svsig = temp("minimap2_pbsv/pbsv/{sample}.svsig.gz"),
    # log:
        # "logs/minimap2_pbsv/pbsv/{sample}.log"
    # shell:
        # """
        # pbsv discover {input.bam} {output.svsig} && \
        # pbsv call {input.genome} {output.svsig} {output.vcf}
        #"""
