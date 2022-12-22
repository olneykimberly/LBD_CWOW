# LBD_CWOW
LBD CWOW RNAseq and DNAseq analysis



### DNA variants
AllSamples.variants.vcf.gz obtained from Dr. Owen Ross's lab.
Subset AllSamples.variants.vcf.gz to only include individuals that are part of this LBD CWOW cohort. 

```
bcftools view -s RES01276,RES01288 AllSamples.variants.vcf.gz > two.sample.test.vcf
```
