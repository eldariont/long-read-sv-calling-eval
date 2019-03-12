import os
import gzip
import sys

include: "rules/mosdepth.smk"
include: "rules/plots.smk"
include: "rules/align.smk"
include: "rules/survivor.smk"
include: "rules/callers.smk"
include: "rules/eval.smk"

configfile: "config.yaml"

##### Target rules #####

rule minimap2:
    input:
        #Alignments
        "minimap2/alignment_stats/alignment_stats.txt",
        "minimap2/mosdepth/regions.combined.gz",
        "minimap2/mosdepth_global_plot/global.html",
        #SV lengths
        expand("minimap2/SV-plots/pooled/SV-length_sniffles_{minsupport}.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("minimap2/SV-plots/pooled/SV-length_svim_default_{max_distance}_{minscore}.png",
                max_distance=[0.1, 0.2, 0.3, 0.4, 0.5], minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        #"minimap2/SV-plots/SV-length_nanosv_pooled.png",
        #Carriers
        expand("minimap2/SV-plots/pooled/SV-sniffles_{minsupport}_carriers.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        #Evaluation
        expand("minimap2/eval/pooled/{run_name}_{max_distance}/tools_pr_svim_sniffles.png", run_name=["default"], max_distance=[0.1, 0.2, 0.3, 0.4, 0.5]),
        expand("minimap2/eval/pooled.subsampled.{fraction}/{run_name}_{max_distance}/tools_pr_svim_sniffles.png", run_name=["default"], max_distance=[0.1, 0.2, 0.3, 0.4, 0.5], fraction=range(10, 91, 10)),
        expand("minimap2/eval/pooled/{run_name}_{max_distance}/svim_pr_multiple_coverages.png", run_name=["default"], max_distance=[0.1, 0.2, 0.3, 0.4, 0.5]),
        expand("minimap2/eval/pooled/{caller}_pr_multiple_coverages.png", caller=["sniffles"])

rule minimap2_pbsv:
    input:
        #Alignments
        "minimap2_pbsv/alignment_stats/alignment_stats.txt",
        "minimap2_pbsv/mosdepth/regions.combined.gz",
        "minimap2_pbsv/mosdepth_global_plot/global.html",
        #SV lengths
        expand("minimap2_pbsv/SV-plots/SV-length_sniffles_{minsupport}_pooled.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("minimap2_pbsv/SV-plots/SV-length_svim_{minscore}_pooled.png",
                minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("minimap2_pbsv/SV-plots/SV-length_pbsv_{minscore}_pooled.png",
                minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"])),
        #Carriers
        expand("minimap2_pbsv/SV-plots/SV-sniffles_{minsupport}_carriers.png",
               minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("minimap2_pbsv/SV-plots/SV-pbsv_{minsupport}_carriers.png",
               minsupport=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"])),
        #Evaluation
        "minimap2_pbsv/eval/pooled/tools_pr_svim_pbsv.png",
        expand("minimap2_pbsv/eval/pooled.subsampled.{fraction}/tools_pr_svim_pbsv.png", fraction=range(10, 91, 10)),
        expand("minimap2_pbsv/eval/pooled/{caller}_pr_multiple_coverages.png", caller=["svim", "pbsv"])

rule minimap2_last_like:
    input:
        expand("minimap2_last_like/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv", "svim"]),
        expand(
            "minimap2_last_like/SV-plots/SV-{caller}_carriers.png",
            caller=["sniffles", "nanosv"]),
        expand("minimap2_last_like/{caller}_combined/sorted_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim"]),
        expand("minimap2_last_like/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        # expand("minimap2_last_like/npinv/{sample}.vcf",
        #       sample=config["samples"]),
        "minimap2_last_like/all_combined/sorted_genotypes.vcf",
        "minimap2_last_like/mosdepth/regions.combined.gz",
        "minimap2_last_like/mosdepth_global_plot/global.html",


rule ngmlr:
    input:
        #Alignments
        "ngmlr/alignment_stats/alignment_stats.txt",
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth_global_plot/global.html",
        #SV lengths
        expand("ngmlr/SV-plots/SV-length_sniffles_{minsupport}_pooled.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("ngmlr/SV-plots/SV-length_svim_{minscore}_pooled.png",
                minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        #"ngmlr/SV-plots/SV-length_nanosv_pooled.png",
        #Carriers
        expand("ngmlr/SV-plots/SV-sniffles_{minsupport}_carriers.png",
               minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        #Evaluation
        "ngmlr/eval/pooled/tools_pr_svim_sniffles.png",
        expand("ngmlr/eval/pooled.subsampled.{fraction}/tools_pr_svim_sniffles.png", fraction=range(10, 91, 10)),
        expand("ngmlr/eval/pooled/{caller}_pr_multiple_coverages.png", caller=["sniffles", "svim"])
