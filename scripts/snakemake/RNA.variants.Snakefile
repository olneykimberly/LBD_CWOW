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
		expand(config["RNA_variants"]+"{sample}_raw_XX.g.vcf.gz", sample = config["female_names"]),

        	expand(config["RNA_variants"]+"{sample}_splitNCigar_XY.bam", sample = config["male_names"]), 	
        	expand(config["RNA_variants"]+"{sample}_XY_recal_data.table", sample = config["male_names"]),	
        	expand(config["RNA_variants"]+"{sample}_bqsr_recal_XY.bam", sample = config["male_names"]), 
        	expand(config["RNA_variants"]+"{sample}_raw_XY.vcf", sample = config["male_names"]),
         	expand(config["RNA_variants"]+"{sample}_raw_snps_XY.vcf", sample = config["male_names"]),
		expand(config["RNA_variants"]+"{sample}_raw_XY.g.vcf.gz", sample = config["male_names"])

#---------------------
# Inferring  variants from RNAseq data 
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
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XX.bam")
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
        recal_table = (config["RNA_variants"]+"{sample}_XX_recal_data.table")
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),
        known_sites = (config["dbsnp138.vcf"]+".gz"),
        gatk = gatk_path
    shell:
        "{params.gatk} BaseRecalibrator -I {input.BAM} -R {params.XX_reference} --known-sites {params.known_sites} -O {output.recal_table}"

# KEY
# BaseRecalibrator detects systematic errors in the estimation of base call accuracy carried out by the sequencing machine
# -I input bam with SplitNCigarReads (see step above)
# -R reference genome used during alignment 
# --known-sites
# -O output table of the recalibration model to be used in the next step (see apply base recall)

# Notes on BaseRecalibrator
# Base quality scores are per-base estimates of error emitted by the sequencing machines; they express how confident the machine was that it called the correct base each time.
# Unfortunately the scores produced by the machines are subject to various sources of systematic (non-random) technical error, leading to over- or under-estimated base quality scores in the data. 
# Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment.
# Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. 
# more information can be found here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-
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
# ApplyBQSR goes through all the reads , using the recalibration file to adjust each base's score based on which bins it falls in
# -I input bam with SplitNCigarReads (see step above)
# -R reference genome used during alignment 
# --bqsr-recal-file the recal_table obtained from BaseRecalibrator step 
# -O output bam file with the base quality scores recalibrated 

# Notes on ApplyBQSR
# new quality score is:
# the sum of the global difference between reported quality scores and the empirical quality
# plus the quality bin specific shift
# plus the cycle x qual and dinucleotide x qual effect
# Following recalibration, the read quality scores are much closer to their empirical scores than before. This means they can be used in a statistically robust manner for downstream processing, such as variant calling.
#---------------------

# we are performing variant calls on germline RNASeq data, thus the “HaplotypeCaller” argument is used.
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
# --java-options -Xmx4g
# HaplotypeCaller variant detection 
# -R reference genome used during alignment 
# -I input bam with recalibrated base quality scores applied (see step above)
# -O output variant call file (vcf)
# -bamout output bam file. Realignment of each haplotype is performed against the reference haplotype to determine potentially variant sites. 

# Notes on HaplotypeCaller
# The program first determines the active regions based on the presence of evidence for variation. 
# Next, for each active region, it builds a De Bruijn-like graph to reassemble the region and identify the haplotypes (a group of alleles in an organism that is inherited together from a single parent) in the data.
# Realignment of each haplotype is performed against the reference haplotype to determine potentially variant sites. 
# Further, for each active region, a pairwise alignment of each read against each haplotype is performed using the PairHMM algorithm – to produce a matrix of likelihoods of haplotypes. 
# These likelihoods are then marginalized to acquire likelihoods of alleles for each potential variant site given the read data. 
# In the end, Bayes’ rule is applied for each potential variant site using the likelihoods of alleles to calculate the likelihoods of each genotype per sample given the read data observed for that sample. 
# eventually, the most likely genotype is assigned to the sample
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

# KEY 
# SelectVariants selecting subsets of variants from a larger variant callset
# -R reference genome used during alignment 
# -V input variant call file 
# --select-type-to-include SNP indiciating the type of variants to select is SNPs 
# -O the output vcf will contain only SNPs 
#--------------------
rule haplotype_gXX:                                                                              
    input:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XX.bam")
    output:                                   
        VCF = (config["RNA_variants"]+"{sample}_raw_XX.g.vcf.gz")                                       
    params:
        XX_reference = (config["GRCh38.Ymasked.fa"]),                                              
        gatk = gatk_path                                                                           
    shell:
        "{params.gatk} --java-options -Xmx4g HaplotypeCaller -R {params.XX_reference} -I {input.BAM} -O {output.VCF} -ERC GVCF"

#---------------------
# XY male samples: 
#---------------------
# format bam files 
rule split_cigar_XY:
    input:
        BAM = (config["starAligned_SCC"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam")
    output:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XY.bam")
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),
        gatk = gatk_path
    shell:
        "{params.gatk} SplitNCigarReads -R {params.XY_reference} -I {input.BAM} -O {output.BAM} || true"

rule base_recall_XY:
    input:
        BAM = (config["RNA_variants"]+"{sample}_splitNCigar_XY.bam")
    output:
        recal_table = (config["RNA_variants"]+"{sample}_XY_recal_data.table")
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

# gvcf
rule haplotype_gXY:                                                                                 
    input:
        BAM = (config["RNA_variants"]+"{sample}_bqsr_recal_XY.bam")                                
    output:                                  
        VCF = (config["RNA_variants"]+"{sample}_raw_XY.g.vcf.gz")                                       
    params:
        XY_reference = (config["GRCh38.YPARs_masked.fa"]),                                         
        gatk = gatk_path                                                                           
    shell:
        "{params.gatk} --java-options -Xmx4g HaplotypeCaller -R {params.XY_reference} -I {input.BAM} -O {output.VCF} -ERC GVCF"
#-----------------------
