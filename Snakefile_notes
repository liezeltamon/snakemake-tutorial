# Comments are advanced options as discussed in main Snakemake tutorial

#SAMPLES = ["A", "B"] # Added to configfile 
configfile: "config.yaml"
 
rule all:
    input:
        "plots/quals.svg"

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map:
    input:
        "data/genome.fa",
        #"data/samples/{sample}.fastq"
        # Using input fuctions:  
        get_bwa_map_input_fastqs
    output:
        #"mapped_reads/{sample}.bam"
	# Specifying output files as tempory i.e. will be deleted if not anymore needed by any rule as input
        temp("mapped_reads/{sample}.bam")
    # Using rule parameters
    params:
	rg=r"@RG\tID:{sample}\tSM:{sample}"
    # Using logging, can also invoke Snakemake with --summary for printing table associating output with rule used to generate it
    log:
        "logs/bwa_mem/{sample}.log"
    # Benchmark directive to measure wall clock time
    benchmark: 
	"benchmarks/{sample}.bwa.benchmark.txt"
    # Specifying threads, if --cores < threads, --cores will be applied
    threads: 8 
    shell:
        #"bwa mem {input} | samtools view -Sb - > {output}"
        # Using rule parameters + threads:
        "bwa mem -R '{params.rg}' -t {threads} {input} | " 
        "samtools view -Sb - > {output} 2> {log}" 
    # Using wrappers from Snakemake wrapper repository
    #wrapper:
    #   "0.15.3/bio/bwa/mem"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        #"sorted_reads/{sample}.bam"
        # Protect important outputs to avoid accidental deletion or modification
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    # Specify environment per rule, then run snakemake --use-conda --cores 1
    #conda: 
	#"envs/samtools.yml"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        # Using config.yml
        #bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),    
        #bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES) 
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf" 
    params:
	rate=config["prior_mutation_rate"]
    logs:
	"logs/bcftools_call/all.log"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
	
