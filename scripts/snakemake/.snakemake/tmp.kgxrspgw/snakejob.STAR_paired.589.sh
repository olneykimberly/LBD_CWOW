#!/bin/sh
# properties = {"type": "single", "rule": "STAR_paired", "local": false, "input": ["../../trimmedReads/NA18-281_FCH5MYFDMXY_L1_trimmed_R1.fastq.gz", "../../trimmedReads/NA18-281_FCH5MYFDMXY_L1_trimmed_R2.fastq.gz"], "output": ["../../starAligned/NA18-281_FCH5MYFDMXY_L1_STAR.bam"], "wildcards": {"sample": "NA18-281_FCH5MYFDMXY_L1"}, "params": {"star": "STAR", "STAR_Index": "/research/labs/neurology/fryer/projects/references/human/GRCh38_def_STAR/", "STAR_GTF": "/research/labs/neurology/fryer/projects/references/human/gencode.v38.annotation.gtf"}, "log": [], "threads": 1, "resources": {}, "jobid": 589, "cluster": {}}
 cd /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake && \
PATH='/research/labs/neurology/fryer/m239830/tools/miniconda3/envs/pigs/bin':$PATH /research/labs/neurology/fryer/m239830/tools/miniconda3/envs/pigs/bin/python3.9 \
-m snakemake ../../starAligned/NA18-281_FCH5MYFDMXY_L1_STAR.bam --snakefile /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw ../../trimmedReads/NA18-281_FCH5MYFDMXY_L1_trimmed_R1.fastq.gz ../../trimmedReads/NA18-281_FCH5MYFDMXY_L1_trimmed_R2.fastq.gz --latency-wait 15 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules STAR_paired --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw/589.jobfinished || (touch /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw/589.jobfailed; exit 1)

