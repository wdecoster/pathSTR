import os
import pandas as pd


work_dir = "/home/wdecoster/pathSTR-1000G/"


# setting up sample dataframes

# VIENNA
samples_vienna = pd.read_table(
    os.path.join(work_dir, "data/all_cram_names_ftp.txt"),
    header=None,
    names=["filename"],
)
samples_vienna_ftp = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/"
)
samples_vienna["sample"] = samples_vienna["filename"].str.replace(".hg38.cram", "")
samples_vienna["path"] = samples_vienna_ftp + samples_vienna["filename"].astype(str)

# MILLER

# SAMPLEINFO
sample_info = pd.read_table(
    os.path.join(work_dir, "data/igsr_samples.tsv"),
    usecols=["Sample name", "Superpopulation code", "Sex"],
).set_index("Sample name")
# Samples with missing Sex or Superpopulation code are dropped
samples = pd.concat([samples_vienna]).set_index("sample").join(sample_info).dropna()

outdir = "/home/wdecoster/pathSTR-1000G"
ref = "/home/wdecoster/database/1KG_ONT_VIENNA_hg38.fa"


def get_path(wildcards):
    return samples.loc[wildcards.sample, "path"]


def get_haploid_chroms(wildcards):
    if samples.loc[wildcards.sample, "Sex"] == "female":
        # all chromosomes are diploid, argument should not be used
        return ""
    else:
        return "--haploid chrX,chrY"


rule all:
    input:
        strdust=expand(
            os.path.join(outdir, "pathSTR_STRdust/{sample}.vcf.gz"),
            sample=samples.index,
        ),
        length_vs_yield=os.path.join(outdir, "plots/length_vs_yield.html"),
        good_samples=os.path.join(
            outdir, "pathSTR_STRdust_good_samples/good_samples.txt"
        ),


rule strdust_unphased:
    output:
        os.path.join(outdir, "pathSTR_STRdust/{sample}.vcf.gz"),
    log:
        "logs/pathSTR_STRdust/{sample}.log",
    params:
        cram=get_path,
        ref=ref,
        targets="/home/wdecoster/p200/data/hg38.STRchive-disease-loci.TRGT.bed",
        binary="/home/wdecoster/repositories/STRdust/target/release/STRdust",
        haploid_chroms=get_haploid_chroms,
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        """RUST_LOG=debug {params.binary} \
        -R {params.targets} \
        --support 2 \
        --unphased \
        --find-outliers \
        {params.haploid_chroms} \
        --somatic \
        {params.ref} \
        {params.cram} 2> {log} \
        | bgzip > {output} 2>> {log}
        """


rule cramino:
    output:
        os.path.join(outdir, "cramino/{sample}.cramino"),
    log:
        "logs/cramino/{sample}.log",
    params:
        cram=get_path,
        ref=ref,
        binary="/home/wdecoster/repositories/cramino/target/release/cramino",
    shell:
        "{params.binary} --karyotype --reference {params.ref} {params.cram} > {output} 2> {log}"


rule cramino_gather:
    input:
        expand(os.path.join(outdir, "cramino/{sample}.cramino"), sample=samples.index),
    output:
        os.path.join(outdir, "cramino/cramino_all.tsv"),
    log:
        "logs/cramino_gather.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(outdir, "scripts/cramino_gather.py"),
    shell:
        "/home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -i {input} -o {output} 2> {log}"


rule plot_length_vs_yield:
    input:
        os.path.join(outdir, "cramino/cramino_all.tsv"),
    output:
        os.path.join(outdir, "plots/length_vs_yield.html"),
    log:
        "logs/plot_length_vs_yield.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(outdir, "scripts/yield_vs_length.py"),
    shell:
        "/home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -i {input} -o {output} 2> {log}"


rule copy_good_samples:
    input:
        overview=os.path.join(outdir, "cramino/cramino_all.tsv"),
        vcfs=expand(
            os.path.join(outdir, "pathSTR_STRdust/{sample}.vcf.gz"),
            sample=samples.index,
        ),
    output:
        os.path.join(outdir, "pathSTR_STRdust_good_samples/good_samples.txt"),
    log:
        "logs/copy_good_samples.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(outdir, "scripts/copy_good_samples.py"),
    shell:
        "outdir=$(dirname {output}) ; /home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -c {input.overview} -v {input.vcfs} -o $outdir 2> {log}"
