# VCF_filter
A multi-threading versatile filter for VCF file, adapt especially for population genetic analysis.

## Usage
``` perl
    --i *input vcf file.
    --o *output file.
    --P *PopMap file.
    --c  coverage for each population.
    --C  total genotyping rate.
    --m  min depth of each individual.
    --M  max depth of each individual.
    --H  max Ho of each population.
    --F  abs Fis value.
    --G  min Genotype quality.
    --Q  min quality for each SNP.
    --g  Global MAF.
    --l  Local MAF.
    --f  apply the filter (Ho and Fis) on each site instead of each population.
    --t  number of threads.
    --s  sort the final vcf file.
    --h  help.
```
