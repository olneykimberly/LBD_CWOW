#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all, kidney, brain and blood sample names
allSamples = list()
type_LBD = list()
LDB_type = list()
sex = list()
read = ["R1", "R2"]
numSamples = 0

with open('sample_read_group_info.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_') # NA00-190.FCH5MYFDMXY_L1_R1_ITCTCTACT-CGCGGTTC.fastq.gz 
        # create a shorter sample name
        stemName = sampleAttributes[0]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "rawReads" : "/research/labs/neurology/fryer/projects/LBD_CWOW/bulkRNA/",
    "rawQC" : "../../rawQC/",
    "trimmedReads" : "../../trimmedReads/",
    "trimmedQC" : "../../trimmedQC/",
    "starAligned" : "../../starAligned/",
    "starAligned_SCC" : "../../startAligned_SCC/",
    "bamstats" : "../../bamstats/",
    "featureCounts" : "../../featureCounts/",
    "kallisto" : "../../kallisto/",
    "salmon" : "../../salmon/", 

    "Comment_Reference" : "This section specifies the location of the Sus scrofa, Ensembl reference genome",
    "GRCh38.fa" : "/research/labs/neurology/fryer/projects/references/human/GRCh38.primary_assembly.genome",
    "GRCh38.star" : "/research/labs/neurology/fryer/projects/references/human/GRCh38_def",
    "GRCh38.gtf" : "/research/labs/neurology/fryer/projects/references/human/gencode.v38.annotation",
    
    "GRCh38.Ymasked.fa" : "/research/labs/neurology/fryer/projects/references/human/GRCh38_Ymasked_XX",
    "GRCh38.YPARs_masked.fa" : "/research/labs/neurology/fryer/projects/references/human/GRCh38_YPARsmasked_XY",

    "GRCh38.Ymasked.star" : "/research/labs/neurology/fryer/projects/references/human/GRCh38_Ymasked_XX",
    "GRCh38.YPARs_masked.star" : "/research/labs/neurology/fryer/projects/references/human/GRCh38_YPARsmasked_XY",

    "dbsnp138.vcf" : "/research/labs/neurology/fryer/projects/references/human/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
    "38.known_indels.vcf" : "/research/labs/neurology/fryer/projects/references/human/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
    "male_names": [],
    "female_names": [],
    "read": {1},
'''
outfile.write(header.format(allSamples, read))

# config formatting
counter = 0
with open('sample_read_group_info.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        sampleName3 = sampleName1.replace("L1","L2")
        sampleName4 = sampleName3.replace("R1","R2")
        base = base.replace("_R1_", "")
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # NA00-190.FCH5MYFDMXY_L1_R1_ITCTCTACT-CGCGGTTC.fastq.gz
        # uniqueNum-number_sequencer_lane_read.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0]
        sampleID = sampleAttributes[0] 
        shortName = stemName
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = sampleID  # Sample
        PU = flowcellID  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/research/labs/neurology/fryer/projects/LBD_CWOW/bulkRNA/",
        "fq1": "{1}",
        "fq2": "{2}",
        "fq3": "{3}",
        "fq4": "{4}",
        "shortName": "{5}",
        "ID": "{6}",
        "SM": "{7}",
        "PU": "{8}",
        "LB": "{9}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, sampleName3, sampleName4, shortName, stemName, sampleID, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()

