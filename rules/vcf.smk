# rule bcftools_reheader_sniffles:
    # """Rule to be deleted as soon as ngmlr uses read groups correctly"""
    # input:
        # "{aligner}/sniffles_genotypes_temp/{sample}.vcf"
    # output:
        # vcf = temp("{aligner}/sniffles_genotypes/{sample}.vcf"),
        # sample = temp("{aligner}/sniffles_genotypes/sample_{sample}.txt")
    # params:
        # sample = "{sample}"
    # log:
        # "logs/{aligner}/bcftools_reheader/{sample}.log"
    # shell:
        # """
        # echo {params.sample} > {output.sample} &&
        # bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
        #"""
