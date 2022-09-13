import os

configfile: "config.json"

#Tools
kallisto_path = "kallisto"

rule all:
        input:
        	expand(config["kallisto"]+"{sample}", sample = config["sample_names"])  


#---------------------
#rule kallisto alignment and quantification: 
#---------------------
rule kallisto_quant:
	input:
		fq1trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		counts = directory((config["kallisto"]+"{sample}"))
	params:
		kallisto = kallisto_path,
		GTF = (config["GRCh38.gtf"]+".gtf"),
		kallisto_index = (config["kallisto_ref_index"]+".fa"),
	shell:
		"{params.kallisto} quant --bias -b 25 -t 8 -i {params.kallisto_index} -g {params.GTF} -o {output.counts} {input.fq1trim} {input.fq2trim} --pseudobam"

# KEY
# quant function to run the quantification algorithm
# --bias learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.
# -b bootstrap 
# -t (nthreads). Specify the number of threads/CPUs used for mapping. (nthreads) The value should be between 1 and 32. 1 by default.
# -i specifics the location and name of the reference index to be used for quantification.
# -g location and name of the reference gene annotation file for transcriptome information. 
# -o name of the directory to write the output to. Unique for each sample ID. 
# pair 1 and pair 2 of the fastq files 
#---------------------
