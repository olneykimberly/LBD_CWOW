import os

configfile: "subset_SCC.RNA.config.json"

#Tools
star_path = "STAR"
picard_path = "picard"
bamtools_path = "bamtools"
gatk_path = "/research/labs/neurology/fryer/m239830/tools/gatk-4.2.6.1/gatk"

rule all:
        input:     
            #expand(config["starAligned_SCC"]+"{sample}_STAR_sort_XX.bam", sample = config["female_names"]),   	   	        		
        	#expand(config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_XX.bam", sample = config["female_names"]),   	
        	expand(config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam", sample = config["female_names"]),   	
        	expand(config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam.bai", sample = config["female_names"]),
        	
        	expand(config["RNA_variants"]+"{sample}_splitNCigar_XX.bam", sample = config["female_names"]), 	
        	expand(config["RNA_variants"]+"{sample}_XX_recal_data.table", sample = config["female_names"]),	
        	expand(config["RNA_variants"]+"{sample}_bqsr_recal_XX.bam", sample = config["female_names"]), 
        	expand(config["RNA_variants"]+"{sample}_raw_XX.vcf", sample = config["female_names"])
  		

#---------------------
# Inferring APOE variants from RNAseq data 
# 	This pipeline follows the recommendations of the GATK pipeline: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
# 		Tutorial: https://expert.cheekyscientist.com/how-to-do-variant-calling-from-rnaseq-ngs-data/ 
#
# Steps
#   1. RNAseq alignment. RNAseq reads are already aligned via STAR two pass mode
#   2. index
#   3. sorting
#   4. marking duplicates
#   5. addig read groups
#   6. bam index
#   7. handling splicing evens in RNAseq data:  gatk SplitNCigarReads -R /research/labs/neurology/fryer/projects/references/human/GRCh38_Ymasked_XX.fa -I NA14-034_FCH5MYFDMXY_L1_STAR_XX.bam -O splitReadsTest.bam
#   8. Base quality recalibration: 
#
#

#---------------------
# XX female samples: 
#---------------------
rule bam_sort_XX:    
    input:
    	IN_BAM = (config["starAligned_SCC"]+"{sample}_STAR_XX.bam")
    output:
        sort_BAM = temporary(config["starAligned_SCC"]+"{sample}_STAR_sort_XX.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_XX:
    input:
        sort_BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_XX.bam")
    output:
        BAM = temporary(config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_XX.bam"),
        metrics = (config["starAligned_SCC"]+"{sample}.picard_sort_mkdup_metrics_XX.txt")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_XX:
    input:
        mkdup_BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_XX.bam")
    output:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    params:
        id = lambda wildcards: config[wildcards.sample]["ID"],
        sm = lambda wildcards: config[wildcards.sample]["SM"],
        lb = lambda wildcards: config[wildcards.sample]["LB"],
        pu = lambda wildcards: config[wildcards.sample]["PU"],
        pl = lambda wildcards: config[wildcards.sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.mkdup_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_XX:
    input:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    output:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam.bai")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

#- GATK
# format bam files 
rule split_cigar_XX:
    input:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    output:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XX.bam")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SplitNCigarReads -R {params.XX_reference} -I {input.BAM} -O {output.BAM} || true"

rule base_recall_XX:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XX.bam")
    output:
        recal_table = (config["RNA_variants"]+"{sample}_XX_recal_data.table")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        known_sites = (config["dbsnp138.vcf"]+".gz"),
        gatk = gatk_path
    shell:
        "{params.gatk} BaseRecalibrator -I {input.BAM} -R {params.XX_reference} --known-sites {params.known_sites} -O {output.recal_table}"
        
rule apply_base_recall_XX:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XX.bam"),
        recal_table = (config["RNA_variants"]+"{sample}_XX_recal_data.table")
    output:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XX.bam")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} ApplyBQSR -R {params.XX_reference} -I {input.BAM} --bqsr-recal-file {input.recal_table} -O {output.BAM}"

# call variants
rule haplotype_XX:
    input:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XX.bam")
    output:
        BAM = (config["RNA_variants"]+"{sample}_bamout_XX.bam"),
        VCF = (config["RNA_variants"]+"{sample}_raw_XX.vcf")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} --java-options -Xmx4g HaplotypeCaller -R {params.XX_reference} -I {input.bam} -O {output.vcf} -bamout {output.BAM}"

#-----------------------
