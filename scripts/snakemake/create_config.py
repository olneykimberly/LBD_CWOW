#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all, kidney, brain and blood sample names
allSamples = list()
brainSamples = list()
kidneySamples = list()
bloodSamples = list()
read = ["R1", "R2"]
numSamples = 0

with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace("-", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # uniqueNum_pigID_tissue_XX_XX_sequencer_lane_read_001.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2] + '_' + \
                   sampleAttributes[5] + '_' + sampleAttributes[6]

        allSamples.append(stemName)
        if stemName.__contains__("Brain"):
            brainSamples.append(stemName)
        elif line.__contains__("Kidney"):
            kidneySamples.append(stemName)
        elif line.__contains__("BB"):
            bloodSamples.append(stemName)
        elif line.__contains__("FB"):
            bloodSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "rawReads" : "/research/labs/neurology/fryer/projects/sepsis/pig/LPS/bulkRNA/",
    "rawQC" : "../../rawQC/",
    "trimmedReads" : "../../trimmedReads/",
    "trimmedQC" : "../../trimmedQC/",
    "starAligned" : "../../starAligned/",
    "bamstats" : "../../bamstats/",
    "featureCounts" : "../../featureCounts/",
    "kallisto" : "../../kallisto/",
    "salmon" : "../../salmon/", 
    "kraken" : "../../kraken/", 
    "kraken_taxo" : "/usr/local/biotools/kraken2/v2.1.1/databases/taxo.k2d",

    "Comment_Reference" : "This section specifies the location of the Sus scrofa, Ensembl reference genome",
    "Sscrofa.Ymasked.fa" : "/research/labs/neurology/fryer/projects/references/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymasked",
    "Scrofa.cdna.Ymasked.fa" : "/research/labs/neurology/fryer/projects/references/pig/Sus_scrofa.Sscrofa11.1.cdna.all.Ymasked",
    "Sscrofa.fa" : "/research/labs/neurology/fryer/projects/references/pig/Sscrofa11.1.dna.toplevel",
    "Sscrofa.gtf" : "/research/labs/neurology/fryer/projects/references/pig/Sus_scrofa.Sscrofa11.1.103",
    "Sscrofa_names_protein.gtf" : "/research/labs/neurology/fryer/projects/references/pig/Sscrofa_geneNames_proteinCoding",
    "salmon_ref_index" : "/research/labs/neurology/fryer/projects/references/pig/salmon_index_Ymasked",
    "star_ref_index" : "/research/labs/neurology/fryer/projects/references/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymasked",
    "kallisto_ref_index" : "/research/labs/neurology/fryer/projects/references/pig/Sus_scrofa.Sscrofa11.1.cdna.all.Ymasked.kallisto",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
    "brain_names": {1},
    "kidney_names": {2},
    "blood_names": {3},
    "read": {4},
'''
outfile.write(header.format(allSamples, brainSamples, kidneySamples, bloodSamples, read))

# config formatting
counter = 0
with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        base = base.replace("_R1_001", "")
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace("-", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # uniqueNum_pigID_tissue_XX_XX_sequencer_lane_read_001.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2] + '_' + \
                   sampleAttributes[5] + '_' + sampleAttributes[6]
        shortName1 = stemName + '_R1'
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
        SM = stemName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/research/labs/neurology/fryer/projects/sepsis/pig/LPS/bulkRNA/",
        "fq1": "{1}",
        "fq2": "{2}",
        "shortName1": "{3}",
        "shortName2": "{4}",
        "ID": "{5}",
        "SM": "{5}",
        "PU": "{6}",
        "LB": "{7}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, shortName1, shortName2, stemName, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
