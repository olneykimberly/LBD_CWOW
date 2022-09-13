import os

configfile: "config.json"

#Tools
#salmon_path = "/research/labs/neurology/fryer/m239830/tools/salmon-1.6.0_linux_x86_64/bin/salmon"
salmon_path = "salmon"

rule all:
        input:
        	expand(config["salmon"]+"{sample}", sample = config["sample_names"])  



#---------------------
#rule salmon alignment and quantification: 
#---------------------
rule salmon_quant:
	input:
		fq1trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		counts = directory((config["salmon"]+"{sample}"))
	params:
		salmon_tool = salmon_path,
		salmon_index = (config["salmon_ref_index"]),
	shell:
		"{params.salmon_tool} quant -p 8 -i {params.salmon_index} -l A -1 {input.fq1trim} -2 {input.fq2trim} -o {output.counts}  --gcBias --seqBias --useVBOpt --writeMappings=salmon.out  --validateMappings --numBootstraps 25 --writeUnmappedNames"

# -i: specify the location of the index directory
# -l A: auto-detect library type. 
# -r: list of files for sample
# -o: output quantification file name
# -p number of threads
# --gcBias to learn and correct for fragment-level GC biases in the input data
# --seqBias will enable Salmon to learn and correct for sequence-specific biases in the input data
# --useVBOpt: use variational Bayesian EM algorithm rather than the ‘standard EM’ to optimize abundance estimates (more accurate)
# --writeMappings=salmon.out: creates a SAM-like file of all of the mappings. We won’t add this parameter, as it creates large files that will take up too much space in your home directory.
# --validateMappings Enables selective alignment of the sequencing reads when mapping them to the transcriptome. This can improve both the sensitivity and specificity of mapping and, as a result, can improve quantification accuracy.
# --numBootstraps compute bootstrapped abundance estimates. This is done by resampling (with replacement) from the counts assigned to the fragment equivalence classes, and then re-running the optimization procedure, either the EM or VBEM, for each such sample. The values of these different bootstraps allows us to assess technical variance in the main abundance estimates we produce. Such estimates can be useful for downstream (e.g. differential expression) tools that can make use of such uncertainty estimates. This option takes a positive integer that dictates the number of bootstrap samples to compute. The more samples computed, the better the estimates of varaiance, but the more computation (and time) required
# --writeUnmappedNames
#---
# more parameter options

# For replicates: 
# salmon quant -i index -l IU -1 <(gunzip -c lib_1_1.fq.gz lib_2_1.fq.gz) -2 <(gunzip -c lib_1_2.fq.gz lib_2_2.fq.gz) --validateMappings -o out
#---------------------

