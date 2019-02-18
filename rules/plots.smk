rule SV_length_plot:
    input:
        "{aligner}/{caller}_calls/pooled.min_{minimum}.vcf"
    output:
        plot = "{aligner}/SV-plots/SV-length_{caller}_{minimum}_pooled.png",
        counts = "{aligner}/SV-plots/SV-nucleotides_affected_{caller}_{minimum}_pooled.txt",
    log:
        "logs/{aligner}/svplot/svlength_{caller}_minimum.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"


rule SV_plot_carriers:
    input:
        "{aligner}/{caller}_calls/pooled.min_{minimum}.vcf"
    output:
        "{aligner}/SV-plots/SV-{caller}_{minimum}_carriers.png"
    log:
        "logs/{aligner}/svplot/svcarriers_{caller}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} --output {output} 2> {log}"
