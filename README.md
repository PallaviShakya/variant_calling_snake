# Variant Calling Protocol 

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

## 4. Run the bash script in slurm 
```
sbatch variant_calling.sh
```