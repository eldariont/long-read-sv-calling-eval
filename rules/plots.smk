rule SV_length_plot:
    input:
        "{aligner}/{caller}_calls/{sample}.min_{minimum}.vcf"
    output:
        plot = "{aligner}/SV-plots/{sample}/SV-length_{caller,sniffles|pbsv}_{minimum}.png",
        counts = "{aligner}/SV-plots/{sample}/SV-nucleotides_affected_{caller,sniffles|pbsv}_{minimum}.txt",
    log:
        "logs/{aligner}/svplot/{sample}/svlength_{caller}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"

rule SV_length_plot_svim:
    input:
        "{aligner}/svim_calls/{sample}/{run_name}_{max_distance}/min_{minimum}.vcf"
    output:
        plot = "{aligner}/SV-plots/{sample}/SV-length_svim_{run_name}_{max_distance}_{minimum}.png",
        counts = "{aligner}/SV-plots/{sample}/SV-nucleotides_affected_svim_{run_name}_{max_distance}_{minimum}.txt",
    log:
        "logs/{aligner}/svplot/{sample}/svlength_svim_{run_name}_{max_distance}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"

rule SV_plot_carriers:
    input:
        "{aligner}/{caller}_calls/{sample}.min_{minimum}.vcf"
    output:
        "{aligner}/SV-plots/{sample}/SV-{caller}_{minimum}_carriers.png"
    log:
        "logs/{aligner}/svplot/{sample}/svcarriers_{caller}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} --output {output} 2> {log}"
