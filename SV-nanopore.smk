import os
import gzip
import sys
import shutil

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
        "minimap2/mosdepth/mean_coverages.txt",
        "minimap2/mosdepth_global_plot/global.html",
        #SV lengths
        expand("minimap2/SV-plots/pooled/SV-length_sniffles_{minsupport}.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("minimap2/SV-plots/pooled/SV-length_pbsv_{minperc}.png",
                minperc=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"])),
        expand("minimap2/SV-plots/pooled/SV-length_svim_default_{max_distance}_{minscore}.png",
                max_distance=[0.3], minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        #"minimap2/SV-plots/SV-length_nanosv_pooled.png",
        #Carriers
        expand("minimap2/SV-plots/pooled/SV-sniffles_{minsupport}_carriers.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("minimap2/SV-plots/pooled/SV-pbsv_{minsupport}_carriers.png",
               minsupport=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"])),
        #Evaluation - Tool comparison
        expand("minimap2/eval/pooled/{run_name}_{max_distance}/tools_pr_all.{vcf}.png", run_name=["default"], max_distance=[0.3], vcf=["vcf", "vcf_het", "vcf.gt"]),
        expand("minimap2/eval/pooled.subsampled.{fraction}/{run_name}_{max_distance}/tools_pr_all.{vcf}.png", run_name=["default"], max_distance=[0.3], fraction=range(10, 91, 10), vcf=["vcf", "vcf_het", "vcf.gt"]),
        #Evaluation - Coverage comparison
        expand("minimap2/eval/pooled/{run_name}_{max_distance}/svim_pr_multiple_coverages.{vcf}.png", run_name=["default"], max_distance=[0.3], vcf=["vcf", "vcf_het", "vcf.gt"]),
        expand("minimap2/eval/pooled/{caller}_pr_multiple_coverages.{vcf}.png", caller=["sniffles", "pbsv"], vcf=["vcf", "vcf_het", "vcf.gt"])


rule ngmlr:
    input:
        #Alignments
        "ngmlr/alignment_stats/alignment_stats.txt",
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth/mean_coverages.txt",
        "ngmlr/mosdepth_global_plot/global.html",
        #SV lengths
        expand("ngmlr/SV-plots/pooled/SV-length_sniffles_{minsupport}.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("ngmlr/SV-plots/pooled/SV-length_svim_default_{max_distance}_{minscore}.png",
                max_distance=[0.3], minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        #Carriers
        expand("ngmlr/SV-plots/pooled/SV-sniffles_{minsupport}_carriers.png",
                minsupport=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        #Evaluation - Tool comparison
        expand("ngmlr/eval/pooled/{run_name}_{max_distance}/tools_pr_all.{vcf}.png", run_name=["default"], max_distance=[0.3], vcf=["vcf", "vcf_het", "vcf.gt"]),
        #expand("ngmlr/eval/pooled.subsampled.{fraction}/{run_name}_{max_distance}/tools_pr_all.{vcf}.png", run_name=["default"], max_distance=[0.3], fraction=range(10, 91, 10), vcf=["vcf", "vcf_het", "vcf.gt"]),
        #Evaluation - Coverage comparison
        #expand("ngmlr/eval/pooled/{run_name}_{max_distance}/svim_pr_multiple_coverages.{vcf}.png", run_name=["default"], max_distance=[0.3], vcf=["vcf", "vcf_het", "vcf.gt"]),
        #expand("ngmlr/eval/pooled/{caller}_pr_multiple_coverages.{vcf}.png", caller=["sniffles"], vcf=["vcf", "vcf_het", "vcf.gt"])
