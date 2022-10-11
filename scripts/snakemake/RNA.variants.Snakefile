import os

configfile: "RNA.config.json"

#Tools
star_path = "STAR"
picard_path = "picard"
bamtools_path = "bamtools"
gatk_path = "/research/labs/neurology/fryer/m239830/tools/gatk-4.2.6.1/gatk"

rule all:
        input:             	
        	expand(config["RNA_variants"]+"{sample}_splitNCigar_XX.bam", sample = config["female_names"]), 	
        	expand(config["RNA_variants"]+"{sample}_XX_recal_data.table", sample = config["female_names"]),	
        	expand(config["RNA_variants"]+"{sample}_bqsr_recal_XX.bam", sample = config["female_names"]), 
        	expand(config["RNA_variants"]+"{sample}_raw_XX.vcf", sample = config["female_names"]),
         	expand(config["RNA_variants"]+"{sample}_raw_snps_XX.vcf", sample = config["female_names"]),

        	expand(config["RNA_variants"]+"{sample}_splitNCigar_XY.bam", sample = config["male_names"]), 	
        	expand(config["RNA_variants"]+"{sample}_XY_recal_data.table", sample = config["male_names"]),	
        	expand(config["RNA_variants"]+"{sample}_bqsr_recal_XY.bam", sample = config["male_names"]), 
        	expand(config["RNA_variants"]+"{sample}_raw_XY.vcf", sample = config["male_names"]),
         	expand(config["RNA_variants"]+"{sample}_raw_snps_XY.vcf", sample = config["male_names"])
	

#---------------------
# Inferring APOE variants from RNAseq data 
# 	This pipeline follows the recommendations of the GATK pipeline: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
# 		Tutorial: https://expert.cheekyscientist.com/how-to-do-variant-calling-from-rnaseq-ngs-data/ 
#
# Steps
#   1. RNAseq alignment. RNAseq reads are already aligned via STAR two pass mode. See RNA.alignment.Snakefile 
#   2. index
#   3. sorting
#   4. marking duplicates
#   5. addig read groups
#   6. bam index
#   7. handling splicing evens in RNAseq data:  gatk SplitNCigarReads -R /research/labs/neurology/fryer/projects/references/human/GRCh38_Ymasked_XX.fa -I NA14-034_FCH5MYFDMXY_L1_STAR_XX.bam -O splitReadsTest.bam
#   8. Base quality recalibration: 
#
# References for base recalibration were downloaded and formatted prior. Follow the steps outlined in this post: https://expert.cheekyscientist.com/how-to-do-variant-calling-from-rnaseq-ngs-data/ 

#---------------------
# XX female samples: 
#---------------------

# split cigar reads
rule split_cigar_XX:
    input:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    output:
        BAM = temporary(config["RNA_variants"]+"{sample}_splitNCigar_XX.bam")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SplitNCigarReads -R {params.XX_reference} -I {input.BAM} -O {output.BAM} || true"

# KEY
# SplitNCigarReads splits reads with N in the cigar into multiple supplementary alignments and hard clips mismatching overhangs. 
# -R reference genome used during alignment
# -I input bam file that is already sorted, duplicates marked, read groups added, and indexed. See RNA.alignment.Snakefile for processing details. 
# -O output bam file with split reads and by default, mapping qualities reassigned to match DNA conventions.
# || true is to allow for zero exit status so snakemake continues to the other steps. 

# Notes on SplitNCigarReads
# we need to reformat some of the alignments that span introns for HaplotypeCaller downstream. 
# the algorithm identifies all N cigars and creates k+1 new reads (where k is the number of N cigar elements). 
# The first read includes the bases to the left of the first N element and the part to the right of the N (including the Ns) is hard clipped.
#---------------------

# base recalibration 
rule base_recall_XX:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XX.bam")
    output:
        recal_table = temporary(config["RNA_variants"]+"{sample}_XX_recal_data.table")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        known_sites = (config["dbsnp138.vcf"]+".gz"),
        gatk = gatk_path
    shell:
        "{params.gatk} BaseRecalibrator -I {input.BAM} -R {params.XX_reference} --known-sites {params.known_sites} -O {output.recal_table}"

# KEY
#---------------------

# apply base recall 
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

# KEY
#---------------------

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
        "{params.gatk} --java-options -Xmx4g HaplotypeCaller -R {params.XX_reference} -I {input.BAM} -O {output.VCF} -bamout {output.BAM}"

# KEY
#---------------------

# select variants
rule select_snps_XX:
    input:
        VCF = (config["RNA_variants"]+"{sample}_raw_XX.vcf")
    output:
        VCF = (config["RNA_variants"]+"{sample}_raw_snps_XX.vcf")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SelectVariants -R {params.XX_reference} -V {input.VCF} --select-type-to-include SNP -O {output.VCF}"

#---------------------
# XY male samples: 
#---------------------
# format bam files 
rule split_cigar_XY:
    input:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam")
    output:
        BAM = temporary(config["RNA_variants"]+"{sample}_splitNCigar_XY.bam")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SplitNCigarReads -R {params.XY_reference} -I {input.BAM} -O {output.BAM} || true"

rule base_recall_XY:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XY.bam")
    output:
        recal_table = temporary(config["RNA_variants"]+"{sample}_XY_recal_data.table")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        known_sites = (config["dbsnp138.vcf"]+".gz"),
        gatk = gatk_path
    shell:
        "{params.gatk} BaseRecalibrator -I {input.BAM} -R {params.XY_reference} --known-sites {params.known_sites} -O {output.recal_table}"
        
rule apply_base_recall_XY:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XY.bam"),
        recal_table = (config["RNA_variants"]+"{sample}_XY_recal_data.table")
    output:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XY.bam")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} ApplyBQSR -R {params.XY_reference} -I {input.BAM} --bqsr-recal-file {input.recal_table} -O {output.BAM}"

# call variants
rule haplotype_XY:
    input:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XY.bam")
    output:
        BAM = (config["RNA_variants"]+"{sample}_bamout_XY.bam"),
        VCF = (config["RNA_variants"]+"{sample}_raw_XY.vcf")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} --java-options -Xmx4g HaplotypeCaller -R {params.XY_reference} -I {input.BAM} -O {output.VCF} -bamout {output.BAM}"

# select variants
rule select_snps_XY:
    input:
        VCF = (config["RNA_variants"]+"{sample}_raw_XY.vcf")
    output:
        VCF = (config["RNA_variants"]+"{sample}_raw_snps_XY.vcf")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SelectVariants -R {params.XY_reference} -V {input.VCF} --select-type-to-include SNP -O {output.VCF}"

#-----------------------
