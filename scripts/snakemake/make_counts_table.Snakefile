import os

configfile: "config.json"


rule all:
	input:
		expand(config["featureCounts"]+"{sample}_countsOnly.txt", sample = config["sample_names"]),
    		expand(config["featureCounts"]+"LBD_81.counts")

#---------------------
#rule reformat counts: 
#---------------------
rule reformat_counts:
    input:
        counts = (config["featureCounts"]+"{sample}_STAR.txt")
    output:
        counts_only = (config["featureCounts"]+"{sample}_countsOnly.txt")
    shell:
        """
        cut -f7 {input.counts} | sed 1d > {output.counts_only};
        sed -i 's,../../starAligned/,,g' {output.counts_only};
        sed -i 's,_STAR_sort_mkdup_rdgrp.bam,,g' {output.counts_only}
        """
  
#---------------------
#rule make tissue counts file: 
#---------------------
rule brain_tissue_counts:
    input:
        counts_only = expand(config["featureCounts"]+"{sample}_countsOnly.txt", sample = config["sample_names"])
    output:
        tissue_counts = config["featureCounts"]+"LBD_81.counts"
    shell:
        """
        paste {input.counts_only} > {output.tissue_counts}
        """
