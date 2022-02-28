import os

configfile: "config.json"

#Tools
fastqc_path = "fastqc"
bbduksh_path = "bbduk.sh"
multiqc_path = "multiqc"
star_path = "STAR"
picard_path = "picard"
bamtools_path = "bamtools"
featureCounts_path = "featureCounts"
kallisto_path = "kallisto"
salmon_path = "salmon"


rule all:
        input:
        	expand(config["rawQC"]+"{sample}_fq1_fastqc.html", sample = config["sample_names"]),
        	expand(config["rawQC"]+"{sample}_fq2_fastqc.html", sample = config["sample_names"]),
        	
        	expand(config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz", sample = config["sample_names"]),
        	expand(config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz", sample = config["sample_names"]),
        	
        	expand(config["trimmedQC"]+"{sample}_trimmed_fq1_fastqc.html", sample = config["sample_names"]),
        	expand(config["trimmedQC"]+"{sample}_trimmed_fq2_fastqc.html", sample = config["sample_names"])
        	        	
        	#expand(config["starAligned"]+"{sample}_STAR.bam", sample = config["sample_names"]),
        	#expand(config["starAligned"]+"{sample}_STAR_sort.bam", sample = config["sample_names"]),    	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup.bam", sample = config["sample_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam", sample = config["sample_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam.bai", sample = config["sample_names"]),  	
        	#expand(config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats.txt", sample = config["sample_names"]),  

        	#expand(config["featureCounts"]+"{sample}_STAR.txt", sample = config["sample_names"])

#---------------------
# Reference genome and annotation were downloaded prior to running snakemake. 
# Ensembl Sscrofa11.1 fasta
# 	wget http://ftp.ensembl.org/pub/release-103/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
# Ensembl gtf
# 	wget http://ftp.ensembl.org/pub/release-103/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.103.gtf.gz   

# prep reference genome. 
# Hard mask Y chromosome for aligning XX samples. 
# Create index and dictionary. 
#---------------------
#rule fastqc on raw:
#---------------------
rule raw_fastqc:
        input:
                fq1 = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq1"] + ".fastq.gz",
                fq2 = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq2"] + ".fastq.gz"
        output:
                fq1_zip =  (config["rawQC"]+"{sample}_fq1_fastqc.zip"),
                fq1_html = (config["rawQC"]+"{sample}_fq1_fastqc.html"),
                fq2_zip =  (config["rawQC"]+"{sample}_fq2_fastqc.zip"),
                fq2_html = (config["rawQC"]+"{sample}_fq2_fastqc.html")
        params:
                fastqc = fastqc_path,
                fastqc_dir = (config["rawQC"]),
                fq1_prefix = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq1"],
                fq2_prefix = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq2"],
        shell:
                """
                {params.fastqc} {input.fq1};
                {params.fastqc} {input.fq2};
                mv {params.fq1_prefix}_fastqc.html {output.fq1_html};
                mv {params.fq1_prefix}_fastqc.zip {output.fq1_zip};
                mv {params.fq2_prefix}_fastqc.html {output.fq2_html};
                mv {params.fq2_prefix}_fastqc.zip {output.fq2_zip}
                """
# KEY
# Run fastqc analysis on read1 and then on read 2. Move the outputs (html and zip) into a new directory
#---------------------
#rule trim fq:
#---------------------
rule trim_bbduk:
	input:
		fq1_trim = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq1"] + ".fastq.gz",
		fq2_trim = lambda wildcards: config[wildcards.sample]["fq_path"] + config[wildcards.sample]["fq2"] + ".fastq.gz"
	output:
		out_fq1 = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		out_fq2 = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} in1={input.fq1_trim} in2={input.fq2_trim} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=/research/labs/neurology/fryer/projects/references/adapters.fa "
		"ktrim=r k=23 mink=11 hdist=1 tpe tbo"	
		
# KEY
# in1/in2 input paired end fastq files
# out1/out2 output paired end fastq files
# ref where adapter fasta is located
# ktrim=r is for right-trimming (3′ adapters), once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left
# ktrim=l is for left-trimming (5′ adapters)
# k=23 kmer length is 23-mers
# mink=11 will additionally look for shorter 11-mers at end of read
# hdist=1 with a small value of mink, it is useful to independently control the hamming/edit distance
# tpe specifies to trim both reads to the same length
# tbo specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences)
#---------------------
#rule fastqc on trimmed:
#---------------------
rule trim_fastqc:
	input:
		fq1_trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2_trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		fq1_zip =  (config["trimmedQC"]+"{sample}_trimmed_fq1_fastqc.zip"),
		fq1_html = (config["trimmedQC"]+"{sample}_trimmed_fq1_fastqc.html"),
		fq2_zip =  (config["trimmedQC"]+"{sample}_trimmed_fq2_fastqc.zip"),
		fq2_html = (config["trimmedQC"]+"{sample}_trimmed_fq2_fastqc.html")
	params:
		fastqc = fastqc_path,
        fastqc_dir = (config["trimmedQC"]),
		fq1_prefix = (config["trimmedReads"]+"{sample}_trimmed_R1"),
		fq2_prefix = (config["trimmedReads"]+"{sample}_trimmed_R2"),
	shell:
		"""
		{params.fastqc} {input.fq1_trim};
		{params.fastqc} {input.fq2_trim};
		mv {params.fq1_prefix}_fastqc.html {output.fq1_html};
		mv {params.fq1_prefix}_fastqc.zip {output.fq1_zip};
		mv {params.fq2_prefix}_fastqc.html {output.fq2_html};
		mv {params.fq2_prefix}_fastqc.zip {output.fq2_zip}
		"""


