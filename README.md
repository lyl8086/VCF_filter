# VCF_filter
A multi-threading versatile filter for VCF file, adapt especially for population genetic analysis.

## Usage
``` perl
    Options [* required, case sensitive]:
    
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
    --p  p-value for hwe.
    --f  apply the filter (Ho and Fis) on each site instead of each population.
    --T  number of threads.
    --s  sort the final vcf file.
    --h  help.
    --totDP total Depth range for eacn site, [min:max].
    --avgDP average Depth range for each site, [min:max].
    --reNb  retain non-bi-allelic sites.
    --reMo  retain monomorphic sites.
    --reIn  retain indels.
```
