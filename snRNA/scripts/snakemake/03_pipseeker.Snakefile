import os
configfile: "config.json"

pipseeker_path = "/research/labs/neurology/fryer/m239830/tools/PIPseeker/pipseeker"

rule all:
        input:
                expand(config["results"]+"/{sample}_results/", sample = config["sample_names"])

#------------------------------
# alignment
#------------------------------
rule pipseeker_full:
        output:
                counts = (directory(config["results"]+"/{sample}_results/")),
        params:
                pipseeker = pipseeker_path, 
                id = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["ID"],
                star_ref = (config["star_ref"])
        shell:
                """
                {params.pipseeker} full --chemistry v3 --fastq {params.id} --star-index-path {params.star_ref} --output-path {output.counts}
                """
                
# $sample is the sampleID (e.g. sample1)
# --full 
# --fastqs is path to the snRNAseq/scRNAseq fastq files
# --star-index-path is the path to the human genome directory. This was created in a prior step.
# --output-path folder to save the outputs

# Pipseeker example tutorial: https://www.fluentbio.com/wp-content/uploads/2023/03/Getting-Started-with-PIPseeker.pdf 
