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
        expand("{trimmed_dir}/{sample}_R2.trimmed.fastq", trimmed_dir=config["trimmed_dir"], sample=SAMPLES)

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
