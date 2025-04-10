# Variant Calling Pipeline

This pipeline can be used to call variants for multiple samples and combine them as one gvcf file and generate genotype, SNP, MNP and INDEL gvcf files. Compared to vcf, gvcf records all the information from a site, so we get all the info whether is a site is variant or not. More info on this format can be found [here.](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-faqs/What_is_a_GVCF_and_how_is_it_different_from_a_%27regular%27_VCF%3F.md) The snakemake pipeline filters based on both bcftools and GATK. GATK is used to flag the variants with low quality and bcftools is used to remove low quality genotypes. 

## 1. Create conda environment using environment.yml 
```
conda env create -f environment.yaml
```

## 2. Index Reference and create sequence dictionary for GATK
```
samtools faids ref.fa
gatk CreateSequenceDictionary -R ref.fa
```

## 3. Edit Snakefile if necessary. 
Pipeline uses GATK and BCFtools to call and filter variants

## 4. Edit config.yaml file 


## 5. Run the bash script in slurm 
```
sbatch variant_calling.sh
```