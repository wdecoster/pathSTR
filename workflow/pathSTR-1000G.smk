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
