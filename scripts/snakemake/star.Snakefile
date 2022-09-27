import os

configfile: "RNA.config.json"

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
        	expand(config["trimmedQC"]+"{sample}_trimmed_fq2_fastqc.html", sample = config["sample_names"]),
        	        	
        	expand(config["starAligned"]+"{sample}_STAR.bam", sample = config["sample_names"]),
        	expand(config["starAligned"]+"{sample}_STAR_sort.bam", sample = config["sample_names"]),    	
        	expand(config["starAligned"]+"{sample}_STAR_sort_mkdup.bam", sample = config["sample_names"]),   	
        	expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam", sample = config["sample_names"]),   	
        	expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bai", sample = config["sample_names"])  	
        	#expand(config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats.txt", sample = config["sample_names"])  

        	#expand(config["starAligned"]+"{sample}_STAR_XX.bam", sample = config["female_names"]),
        	#expand(config["starAligned"]+"{sample}_STAR_sort_XX.bam", sample = config["female_names"]), 
        	#expand(config["starAligned"]+"{sample}_STAR_sort_XX.bam.bai", sample = config["female_names"]),  	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_XX.bam", sample = config["female_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam", sample = config["female_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam.bai", sample = config["female_names"]),  	
        	#expand(config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats_XX.txt", sample = config["female_names"]),  

        	#expand(config["starAligned"]+"{sample}_STAR_XY.bam", sample = config["male_names"]),
        	#expand(config["starAligned"]+"{sample}_STAR_sort_XY.bam", sample = config["male_names"]),  
        	#expand(config["starAligned"]+"{sample}_STAR_sort_XY.bam.bai", sample = config["male_names"]),  	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_XY.bam", sample = config["male_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam", sample = config["male_names"]),   	
        	#expand(config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam.bai", sample = config["male_names"]),  	
        	#expand(config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats_XY.txt", sample = config["male_names"])  

#---------------------
# Reference genome and annotation were downloaded prior to running snakemake. 
# prep reference genome, create index and dictionary. 

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

#---------------------
#rule star alignment: 
#---------------------

rule STAR_paired:
	input:
		fq1_trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2_trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		out_1 = (config["starAligned"]+"{sample}_STAR.bam")
	params:
		star = star_path,
		STAR_Index = (config["GRCh38.star"]+"_STAR/"),
		STAR_GTF = (config["GRCh38.gtf"]+".gtf"),
	shell:
		"""
		{params.star} --runThreadN 8 --genomeDir {params.STAR_Index} --sjdbGTFfile {params.STAR_GTF} --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat --readFilesIn {input.fq1_trim} {input.fq2_trim} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1};
		mv {output.out_1}Aligned.out.bam {output.out_1}
		"""

# KEY
#--runThreadN NumberOfThreads
#--genomeDir specifies where indices are located
#--sjdbGTFfile gene annotation file, used for splice aware junctions
#--twopassMode TAR will perform the 1st pass mapping, then it will automatically extract junctions, insert them into the genome index, and, finally, re-map all reads in the 2nd mapping pass.
#--readFilesCommand zcat for reading in compressed .gz files 
#--readFilesIn read in pair end trimmed fastq files
#--outSAMtype BAM Unsorted. Output will not be sorted by coordinate
#--quantMode TranscriptomeSAM GeneCounts to get the aligned transcripts an counts the number of reads per gene id
#--outFileNamePrefix. Naming prefix for the output bam file. 

# Output 
# SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that STAR defines the junction start/end as intronic bases, while many other software define them as exonic bases.

# quantMode GeneCounts
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)


# Note that STAR counts a paired-end read as one read, (unlike the samtools flagstat/idxstats, which count each mate separately). 
#--quantMode geneCounts, can be used to get gene counts 
#--genomeLoad LoadAndExit. Loads the star genome index into memory to be used by all star jobs. Will unload after script is done running. 
#--sjdbOverhang specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.
#---------------------
#rule bam sort: 
#---------------------
rule STAR_bam_sort:    
    input:
    	IN_BAM = (config["starAligned"]+"{sample}_STAR.bam")
    output:
        sort_BAM = temporary(config["starAligned"]+"{sample}_STAR_sort.bam")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g SortSam I={input.IN_BAM} O={output.sort_BAM} SORT_ORDER=coordinate"

#---- bamtools sorting alternative
#        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"
# KEY
# sort. Sort bam by coordinates
# -in. input bam file
# -out. output bam file
#---------------------
#rule mark duplicates sort: 
#---------------------
rule STAR_MarkDups:
    input:
        sort_BAM = (config["starAligned"]+"{sample}_STAR_sort.bam")
    output:
        BAM = temporary(config["starAligned"]+"{sample}_STAR_sort_mkdup.bam"),
        metrics = (config["starAligned"]+"{sample}.picard_sort_mkdup_metrics.txt")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

# KEY
#-Xmx14g. -Xms<size> it' initial heap size. 14GB for the heap
#MarkDuplicates. This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.
#I=. input bam file
#O=. output bam file
#M=. File to write duplication metrics to
#VALIDATION_STRINGENCY=LENIENT. Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
#---------------------
#rule add read groups sort: 
#---------------------
rule STAR_AddReadGrps:
    input:
        mkdup_BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam")
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
# KEY
#-Xmx14g. -Xms<size> it' initial heap size. 14GB for the heap
#AddOrReplaceReadGroups
#I=. input bam file
#O=. output bam file
#RGID = Read group identifier 
#GPU = Platform Unit 
#GSM = Sample 
#GPL = Platform/technology used to produce the read 
#GLB = DNA preparation library identifier 
#VALIDATION_STRINGENCY=LENIENT. Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
#---------------------
#rule index bam: 
#---------------------
rule STAR_index_bam:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bai")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g BuildBamIndex -I {input.BAM}"

#------ bamtools index alternative
# "{params.bamtools} index -in {input.BAM}"

# KEY
#index. Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access.
#-in. input bam file

# many downstream tools require the bam to be indexed. 
#---------------------
#rule bam stats: 
#---------------------
rule STAR_stats_bam:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp.bam")
    output:
        stats = (config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats.txt")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

# KEY
#stats. Prints some basic statistics from input BAM file(s)
#-in. input bam file
# stats output will include number of reads mapped and number of duplicate reads.


#-----------------------------------------------------------------------------------------
# Female XX only samples
#-----------------------------------------------------------------------------------------

rule STAR_paired_XX:
	input:
		fq1_trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2_trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		out_1 = (config["starAligned"]+"{sample}_STAR_XX.bam")
	params:
		star = star_path,
		STAR_Index = (config["GRCh38.Ymasked.star"]+"_STAR/"),
		STAR_GTF = (config["GRCh38.gtf"]+".gtf"),
	shell:
		"""
		{params.star} --runThreadN 8 --genomeDir {params.STAR_Index} --sjdbGTFfile {params.STAR_GTF} --readFilesCommand zcat --readFilesIn {input.fq1_trim} {input.fq2_trim} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1};
		mv {output.out_1}Aligned.out.bam {output.out_1}
		"""

rule STAR_bam_sort_XX:    
    input:
    	IN_BAM = (config["starAligned"]+"{sample}_STAR_XX.bam")
    output:
        sort_BAM = temporary(config["starAligned"]+"{sample}_STAR_sort_XX.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"


rule STAR_MarkDups_XX:
    input:
        sort_BAM = (config["starAligned"]+"{sample}_STAR_sort_XX.bam")
    output:
        BAM = temporary(config["starAligned"]+"{sample}_STAR_sort_mkdup_XX.bam"),
        metrics = (config["starAligned"]+"{sample}.picard_sort_mkdup_metrics_XX.txt")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"


rule STAR_AddReadGrps_XX:
    input:
        mkdup_BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_XX.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
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

rule STAR_index_bam_XX:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam.bai")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"


rule STAR_stats_bam_XX:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XX.bam")
    output:
        stats = (config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats_XX.txt")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


#-----------------------------------------------------------------------------------------
# Male XY only samples
#-----------------------------------------------------------------------------------------
rule STAR_paired_XY:
	input:
		fq1_trim = (config["trimmedReads"]+"{sample}_trimmed_R1.fastq.gz"),
		fq2_trim = (config["trimmedReads"]+"{sample}_trimmed_R2.fastq.gz")
	output:
		out_1 = (config["starAligned"]+"{sample}_STAR_XY.bam")
	params:
		star = star_path,
		STAR_Index = (config["GRCh38.YPARs_masked.star"]+"_STAR/"),
		STAR_GTF = (config["GRCh38.gtf"]+".gtf"),
	shell:
		"""
		{params.star} --runThreadN 8 --genomeDir {params.STAR_Index} --sjdbGTFfile {params.STAR_GTF} --readFilesCommand zcat --readFilesIn {input.fq1_trim} {input.fq2_trim} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1};
		mv {output.out_1}Aligned.out.bam {output.out_1}
		"""

rule STAR_bam_sort_XY:    
    input:
    	IN_BAM = (config["starAligned"]+"{sample}_STAR_XY.bam")
    output:
        sort_BAM = temporary(config["starAligned"]+"{sample}_STAR_sort_XY.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"


rule STAR_MarkDups_XY:
    input:
        sort_BAM = (config["starAligned"]+"{sample}_STAR_sort_XY.bam")
    output:
        BAM = temporary(config["starAligned"]+"{sample}_STAR_sort_mkdup_XY.bam"),
        metrics = (config["starAligned"]+"{sample}.picard_sort_mkdup_metrics_XY.txt")
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"


rule STAR_AddReadGrps_XY:
    input:
        mkdup_BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_XY.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam")
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

rule STAR_index_bam_XY:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam")
    output:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam.bai")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"


rule STAR_stats_bam_XY:
    input:
        BAM = (config["starAligned"]+"{sample}_STAR_sort_mkdup_rdgrp_XY.bam")
    output:
        stats = (config["bamstats"]+"{sample}_STAR_sort_mkdup_rdgrp_stats_XY.txt")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
#---------------------
# End of Snakefile. Proceed to R scripts for differential expression. 
