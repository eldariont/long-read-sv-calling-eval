import os
import gzip
import sys
import shutil

configfile: "config/config.yaml"

ALIGNERS=["minimap2"]
SUBSAMPLES=["pooled"] + [("pooled.subsampled." + str(fraction)) for fraction in range(10, 100, 10)]
VCFS=["giab", "giab.gt"]

wildcard_constraints:
    aligner="minimap2|ngmlr|pbmm2"

include: "workflow/rules/mosdepth.smk"
include: "workflow/rules/plots.smk"
include: "workflow/rules/align.smk"
#include: "workflow/rules/survivor.smk"
include: "workflow/rules/callers.smk"
include: "workflow/rules/eval.smk"

##### Target rules #####

rule all:
    input:
        #Alignments
        expand("pipeline/alignment_stats/alignment_stats.{aligner}.txt", aligner=ALIGNERS),
        #"minimap2/mosdepth/regions.combined.gz",
        #"minimap2/mosdepth/mean_coverages.txt",
        #"minimap2/mosdepth_global_plot/global.html",
        #SV lengths
        #"pipeline/SV-plots/minimap2/pooled/SV-length_Sniffles_11.png",
        #Evaluation
        expand("pipeline/eval/{aligner}/results.{aligner}.all.png", aligner=ALIGNERS),
        expand("pipeline/eval/{aligner}/results.{aligner}.tools.{vcf}.png", aligner=ALIGNERS, vcf=VCFS),
        expand("pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.png", aligner=ALIGNERS, vcf=VCFS),
        expand("pipeline/eval/{aligner}/results.{aligner}.svim.parameters.png", aligner=ALIGNERS),
        #expand("pipeline/eval/{aligner}/results.{aligner}.coverages.bar.png", aligner=ALIGNERS),
        expand("pipeline/alignment_stats/alignment_stats.{aligner}.txt", aligner=ALIGNERS),
        expand("pipeline/mosdepth/mean_coverages.{aligner}.txt", aligner=ALIGNERS)
