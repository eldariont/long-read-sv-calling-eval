localrules: SV_length_plot, SV_length_plot_svim

rule SV_length_plot:
    input:
        "pipeline/{caller}/{aligner}/{data}/min_{minimum}.indel.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_{caller,Sniffles|pbsv}_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-nucleotides_affected_{caller,sniffles|pbsv}_{minimum}.txt",
    log:
        "logs/svplot/svlength_{caller}_{aligner}_{data}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"

rule SV_length_plot_svim:
    input:
        "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minimum}.indel.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_SVIM_{parameters}_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-nucleotides_affected_SVIM_{parameters}_{minimum}.txt",
    log:
        "logs/svplot/svlength_SVIM_{aligner}_{data}_{parameters}_{minimum}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"