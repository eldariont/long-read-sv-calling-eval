rule SV_length_plot:
    input:
        "{aligner}/{caller}_calls/pooled.vcf"
    output:
        plot = "{aligner}/SV-plots/SV-length_{caller}_pooled.png",
        counts = "{aligner}/SV-plots/SV-nucleotides_affected_{caller}_pooled.txt",
    log:
        "logs/{aligner}/svplot/svlength_{caller}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"


rule SV_plot_carriers:
    input:
        "{aligner}/{caller}_calls/pooled.vcf"
    output:
        "{aligner}/SV-plots/SV-{caller}_carriers.png"
    log:
        "logs/{aligner}/svplot/svcarriers_{caller}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} --output {output} 2> {log}"
