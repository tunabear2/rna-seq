configfile: "config/config.yaml"

import os
from glob import glob
import yaml
from snakemake.io import expand, glob_wildcards

rowdata_dir = config["rowdata_dir"]

SAMPLES, = glob_wildcards(os.path.join(rowdata_dir, "{sample}_1.fastq"))

def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)
def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)

rule all:
    input:
        expand("{fastqc_dir}/{sample}_1_fastqc.zip", fastqc_dir=config["fastqc_dir"], sample=SAMPLES),
        expand("{fastqc_dir}/{sample}_2_fastqc.zip", fastqc_dir=config["fastqc_dir"], sample=SAMPLES),
        expand("{fastqc_dir}/{sample}_1_fastqc.html", fastqc_dir=config["fastqc_dir"], sample=SAMPLES),
        expand("{fastqc_dir}/{sample}_2_fastqc.html", fastqc_dir=config["fastqc_dir"], sample=SAMPLES),
        expand("{trimmed_dir}/{sample}_R1.trimmed.fastq", trimmed_dir=config["trimmed_dir"], sample=SAMPLES),
        expand("{trimmed_dir}/{sample}_R2.trimmed.fastq", trimmed_dir=config["trimmed_dir"], sample=SAMPLES),
        expand("results/star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/rsem/{sample}.genes.results", sample=SAMPLES),
        expand("results/rsem/{sample}.isoforms.results", sample=SAMPLES)

rule fastqc:
    input:
        r1 = lambda wildcards: f"data/fastq/{wildcards.sample}_1.fastq",
        r2 = lambda wildcards: f"data/fastq/{wildcards.sample}_2.fastq",
        direc = config["fastqc_dir"]
    output:
        r1_qc = "results/fastqc/{sample}_1_fastqc.zip",
        r2_qc = "results/fastqc/{sample}_2_fastqc.zip",
        r1_html = "results/fastqc/{sample}_1_fastqc.html",
        r2_html = "results/fastqc/{sample}_2_fastqc.html"
    singularity:
        config["tools"]["fastqc"]
    shell:
        "fastqc {input.r1} {input.r2} -o {input.direc}"
rule fastp:
    input:
        r1 = lambda wildcards: f"data/fastq/{wildcards.sample}_1.fastq",
        r2 = lambda wildcards: f"data/fastq/{wildcards.sample}_2.fastq"
    output:
        r1_trimmed = "results/fastp/{sample}_R1.trimmed.fastq",
        r2_trimmed = "results/fastp/{sample}_R2.trimmed.fastq",
        html = "results/fastp/{sample}.html",
        json = "results/fastp/{sample}.json"
    singularity:
        config["tools"]["fastp"]
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1_trimmed} -O {output.r2_trimmed} --html {output.html} --json {output.json}"
rule star_index:
    input:
        ref = "data/reference/GRCh38.p14.genome.fa",
        gtf = "data/reference/gencode.v49.annotation.gtf"
    output:
        idx = "data/star_idx/"
    singularity:
        config["tools"]["star"]
    shell:
        "STAR --runThreadN 6 --runMode genomeGenerate --genomeDir {output.idx} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} --sjdbOverhang 99"
rule star:
    input:
        idx = "data/star_idx/",
        r1 = lambda wildcard: f"results/fastp/{wildcard.sample}_R1.trimmed.fastq",
        r2 = lambda wildcard: f"results/fastp/{wildcard.sample}_R2.trimmed.fastq"
    output:
        bam = "results/star/{sample}.Aligned.sortedByCoord.out.bam"
    params:
        prefix = lambda wc: f"results/star/{wc.sample}."
    singularity:
        config["tools"]["star"]
    shell:
        "STAR --runThreadN 6 --genomeDir {input.idx} --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode TranscriptomeSAM --twopassMode None"
rule rsem_ref:
    input:
        ref = "data/reference/GRCh38.p14.genome.fa",
        gtf = "data/reference/gencode.v49.annotation.gtf"
    output:
        directory("data/rsem")
    singularity:
        config["tools"]["rsem"]
    shell:
        "rsem-prepare-reference --gtf {input.gtf} {input.ref} data/rsem/rsem_genv49_GR38p14"
rule rsem:
    input:
        bam_tx = lambda wildcards: f"results/star/{wildcards.sample}.Aligned.toTranscriptome.out.bam"
    output:
        gene = "results/rsem/{sample}.genes.results",
        iso = "results/rsem/{sample}.isoforms.results"
    singularity:
        config["tools"]["rsem"]
    params:
        output_prefix = lambda wc: f"results/rsem/{wc.sample}"
    shell:
        "rsem-calculate-expression --paired-end --alignments --fai data/reference/GRCh38.p14.genome.fa.fai --strandedness none --estimate-rspd --calc-pme --no-bam-output {input.bam_tx} data/rsem/rsem_genv49_GR38p14 {params.output_prefix}"
