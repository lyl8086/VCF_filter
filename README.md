# VCF_filter
A multi-threading versatile filter for VCF file, adapt especially for population genetic analysis.

## Usage
- --i [required]
    input vcf/bcf file, gzipped accepted.
- --o [required]
    output vcf file.
- --P [required]
    a popMap file with two colums, corresponding to individual and population name, separated by tab, see [popmap.txt](samples/popmap.txt).
- --c
    individual coverage for each population, can be a ratio (<=1) or number (>1).
- --C
    total genotyping rate for a site, can be a ratio (<=1) or number (>1).
- --m
    minimum depth for each individual genotype.
- --M
    maximum depth for a individual genotype.
- --H
    maximum observed heterozygosity (Ho) of each population.
- --F
    absolute maximum inbreeding coefficient (Fis) of each population.
- --G
    minimum Genotype quality (GQ).
- --Q
    minimum quality for each SNP (Q).
- --g
    global minor allele frequncy for a site.
- --l
    local minor allele frequncy, if a population has this many minor allele, evaluated by the given value, then retain this site.
- --p
    sites with a p-value below the threshold defined by this option are taken to be out of HWE, and therefore excluded.
- --f
    apply the filter (--H and --F) on each site instead of each population.
- --T
    number of threads.
- --s
    sort the final vcf file, recommond on.
- --h
    help message.
- --totDP
    total Depth range for each site, [min:max], e.g., 10:100 means the total depth for this site should between 10 and 100.
- --avgDP
    average Depth range for each site, [min:max], e.g., 10:100 means the average depth for this site should between 10 and 100.
- --reNb
    retain non-bi-allelic sites, default unset.
- --reMo
    retain monomorphic sites, default unset.
- --reIn  
    retain indels, defalt unset.
- --rePo
    threshold to define a polymorphic loci of each populations, should be a ratio.


## Example usage

```shell
VCF_filter_multi.pl --i samples.vcf.gz --P popmap.txt --o qc.vcf.gz --c 0.8 --C 0.9 --m 6 --H 0.75 --G 20 --Q 20 --g 0.05 --l 0.2 --p 0.01 --avgDP 6:20 --T 8 --s --f

Which set:
at least 80% individuals of each population cover this site
at least 90% individuals totally cover this site
minimum read depth of 6
maximum observed heterozygosity considering all individuals (as --f set) <=0.75
genotyping quality >=20
snp quality >=20
minor allele frequency >=5%
retain a site if one population with minor allele frequency >=20%
exclude a site if p-value of HWE test <=0.01
average read depth for a site should between 6 and 20
use 8 cpus
sort the final vcf
```

