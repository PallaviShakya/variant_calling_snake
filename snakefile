#Loading samples from config file

configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())

rule all:
    input:
        # expand("results_F2/fastqc/{sample}_R1.fastqc.html", sample = SAMPLES),
        # expand("results_F2/fastqc/{sample}_R2.fastqc.html", sample = SAMPLES),
        expand("results_F2/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        "results_F2/combined/combined.g.vcf.gz",
        "results_F2/variants/combined_genotyped.vcf.gz",
        "results_F2/var/combined_SNP.vcf",
        "results_F2/var/combined_INDEL.vcf",
        "results_F2/var/combined_MNP.vcf",
        "results_F2/var/combined_SNP_genotypes.txt",
        "results_F2/variants/combined_genotyped.filtered.vcf.gz",
        "results_F2/variants/combined_genotyped.filtered.clean.vcf.gz",
        "results_F2/var/combined_SNP_filtered.vcf",
        "results_F2/var/combined_INDEL_filtered.vcf",
        "results_F2/var/combined_MNP_filtered.vcf",
        "results_F2/var/combined_SNP_genotypes_filtered.txt"

        # expand("results_F2/var/{sample}_SNP.vcf", sample=SAMPLES),
        # expand("results_F2/var/{sample}_INDEL.vcf", sample=SAMPLES),
        # expand("results_F2/var/{sample}_MNP.vcf", sample= SAMPLES),
        # expand("results_F2/var/{sample}_SNP_genotypes.txt", sample=SAMPLES)


#QC illumina reads

rule fastqc:
    input: 
        R1=lambda wildcards: config['samples'][wildcards.sample]['R1'],
        R2=lambda wildcards: config['samples'][wildcards.sample]['R2']
    output: 
        "results_F2/fastqc/{sample}_R1.fastqc.html",
        "results_F2/fastqc/{sample}_R2.fastqc.html"
    log:
        "logs_for_F2/fastqc_{sample}.log"
    shell:
        "fastqc {input} -o results_F2/fastqc > {log} 2>&1"



# #Mapping Illumina reads to the reference genome 

rule bwa_mem: 
    input: 
        ref=config['reference'],
        R1=lambda wildcards: config['samples'][wildcards.sample]['R1'],
        R2=lambda wildcards: config['samples'][wildcards.sample]['R2']
    output: 
        "results_F2/mapped/{sample}.sorted.bam"
    log:
        "logs_for_F2/bwa_mem_{sample}.log"
        
    shell: 
        """
        bwa mem {input.ref} {input.R1} {input.R2}| \
        samtools view -bS - | \
        samtools sort -o {output} > {log} 2>&1
        """

rule index_bam: 
    input: 
        "results_F2/mapped/{sample}.sorted.bam"
    output: 
        "results_F2/mapped/{sample}.sorted.bam.bai"
    log:
        "logs_for_F2/bam_index.{sample}.sorted.log"
    shell:
        "samtools index {input} > {log} 2>&1"

# # # # Step to add read groups using Picard
rule add_read_groups:
    input:
        bam="results_F2/mapped/{sample}.sorted.bam"
    output:
        bam_rg="results_F2/mapped/{sample}.sorted.rg.bam"
    log:
        "logs_for_F2/add_read_groups_{sample}.log"
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam_rg} \
            RGID=1 \
            RGLB=lib1 \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM={wildcards.sample}> {log} 2>&1
            
        """

rule index_bam_2:
    input:
        "results_F2/mapped/{sample}.sorted.rg.bam"
    output:
        "results_F2/mapped/{sample}.sorted.rg.bam.bai"
    log:
        "logs_for_F2/bam_index_{sample}.log"
    shell:
        "samtools index {input} > {log} 2>&1"


# # # #Variant calling: 

rule variant_calling: 
    input: 
        bam="results_F2/mapped/{sample}.sorted.rg.bam",
        bai="results_F2/mapped/{sample}.sorted.rg.bam.bai",
        ref=config['reference']
    output:
        "results_F2/variants/{sample}.vcf"
    log: 
        "logs_for_F2/variant_calling_{sample}.log"
    shell: 
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -ploidy 2 > {log} 2>&1

        """

rule haplotypecaller_gvcf: 
    input: 
        bam="results_F2/mapped/{sample}.sorted.rg.bam",
        bai="results_F2/mapped/{sample}.sorted.rg.bam.bai",
        ref=config["reference"]
    output:
        "results_F2/gvcf/{sample}.g.vcf.gz"
    log:
        "logs_for_F2/haplotypecaller_{sample}.log"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -ERC GVCF > {log} 2>&1
        """

# # Merge per sample gvcfs per sample
rule combine_gvcfs: 
    input: 
        gvcfs=expand("results_F2/gvcf/{sample}.g.vcf.gz", sample=SAMPLES)
    output:
        "results_F2/combined/combined.g.vcf.gz"
    log:
        "logs_for_F2/combine_gvcfs.log"
    params:
        variants=lambda wildcards, input: " ".join(f"--variant {g}" for g in input.gvcfs)
    shell:
        """
        gatk CombineGVCFs \
            -R {config[reference]} \
            {params.variants} \
            -O {output} > {log} 2>&1

        """
    
# # Joint genotyping
rule genotype_gvcfs:
    input:
        "results_F2/combined/combined.g.vcf.gz"
    output:
        "results_F2/variants/combined_genotyped.vcf.gz"
    log:
        "logs_for_F2/genotypegvcfs.log"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {config[reference]} \
            -V {input} \
            -O {output} > {log} 2>&1
        """



# # Calling SNPs: 

rule snp_call: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.vcf.gz"
    output:
        "results_F2/var/combined_SNP.vcf"
    log:
        "logs_for_F2/snp_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include SNP \
            -O {output} > {log} 2>&1
        """

rule indel_call: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.vcf.gz"
    output:
        "results_F2/var/combined_INDEL.vcf"
    log:
        "logs_for_F2/indel_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include INDEL \
            -O {output} > {log} 2>&1
        """

# #Calling MNP: 

rule mnp_call: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.vcf.gz"
    output:
        "results_F2/var/combined_MNP.vcf"
    log:
        "logs_for_F2/mnp_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include MNP \
            -O {output} > {log} 2>&1
        """

# #extracting genotypes from SNP
rule call_geno: 
    input:
        "results_F2/var/combined_SNP.vcf"
    output:
        "results_F2/var/combined_SNP_genotypes.txt"
    log:
        "logs_for_F2/extract_geno_combined.log"
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F GT \
            -O {output} > {log} 2>&1
        """

# Variant with filters to flag the low quality variants
rule variant_filtering:
    input:
        vcf="results_F2/variants/combined_genotyped.vcf.gz",
        ref=config["reference"]
    output:
        "results_F2/variants/combined_genotyped.filtered.vcf.gz"
    log:
        "logs_for_F2/variant_filtering.log"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --filter-name "LOW_QUAL" \
            -O {output} > {log} 2>&1
        """

# bcf tools to clean according to GQ and DP

rule filter_genotypes:
    input:
        "results_F2/variants/combined_genotyped.filtered.vcf.gz"
    output:
        "results_F2/variants/combined_genotyped.filtered.clean.vcf.gz"
    log:
        "logs_for_F2/genotype_filtering.log"
    shell:
        """
        bcftools view \
            -f PASS \
            -i 'FMT/GQ > 20 & FMT/DP > 10' \
            {input} -Oz -o {output} > {log} 2>&1
        tabix -p vcf {output}
        """


rule snp_call_after_filter: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.filtered.clean.vcf.gz"
    output:
        "results_F2/var/combined_SNP_filtered.vcf"
    log:
        "logs_for_F2/snp_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include SNP \
            --exclude-filtered \
            -O {output} > {log} 2>&1
        """

rule indel_call_filtered: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.filtered.clean.vcf.gz"
    output:
        "results_F2/var/combined_INDEL_filtered.vcf"
    log:
        "logs_for_F2/indel_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include INDEL \
            --exclude-filtered \
            -O {output} > {log} 2>&1
        """

# #Calling MNP: 

rule mnp_call_filtered: 
    input:
        ref=config['reference'],
        var="results_F2/variants/combined_genotyped.filtered.clean.vcf.gz"
    output:
        "results_F2/var/combined_MNP_filtered.vcf"
    log:
        "logs_for_F2/mnp_calling_combined.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.var} \
            --select-type-to-include MNP \
            --exclude-filtered \
            -O {output} > {log} 2>&1
        """

# #extracting genotypes from SNP
rule call_geno_filtered: 
    input:
         "results_F2/var/combined_SNP_filtered.vcf"
    output:
        "results_F2/var/combined_SNP_genotypes_filtered.txt"
    log:
        "logs_for_F2/extract_geno_combined_filtered.log"
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F GT \
            -O {output} > {log} 2>&1
        """