#!/bin/sh
# properties = {"type": "single", "rule": "STAR_bam_sort", "local": false, "input": ["../../starAligned/NA13-047_FCH5MYFDMXY_L2_STAR.bam"], "output": ["../../starAligned/NA13-047_FCH5MYFDMXY_L2_STAR_sort.bam"], "wildcards": {"sample": "NA13-047_FCH5MYFDMXY_L2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 712, "cluster": {}}
 cd /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake && \
PATH='/research/labs/neurology/fryer/m239830/tools/miniconda3/envs/pigs/bin':$PATH /research/labs/neurology/fryer/m239830/tools/miniconda3/envs/pigs/bin/python3.9 \
-m snakemake ../../starAligned/NA13-047_FCH5MYFDMXY_L2_STAR_sort.bam --snakefile /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw ../../starAligned/NA13-047_FCH5MYFDMXY_L2_STAR.bam --latency-wait 15 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules STAR_bam_sort --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw/712.jobfinished || (touch /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/.snakemake/tmp.kgxrspgw/712.jobfailed; exit 1)

