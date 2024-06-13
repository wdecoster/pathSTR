import os
import pandas as pd


def parse_motifs(info):
    """Return the length of the longest motif"""
    motif_string = [i for i in info.split(";") if i.startswith("MOTIFS=")][0]
    motifs = motif_string.split("=")[1].split(",")
    return max([len(m) for m in motifs])


# setting up paths
work_dir = "/home/wdecoster/pathSTR-1000G/"
strdust = "/home/wdecoster/repositories/STRdust/target/release/STRdust"
longTR = "/home/wdecoster/repositories/LongTR/LongTR"  # hacked version of LongTR to enable remote cram
strchive = pd.read_csv(
    "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/hg38.STRchive-disease-loci.TRGT.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "info"],
)
targets_strdust = os.path.join(work_dir, "data/STRchive_STRdust.bed")
strchive.to_csv(targets_strdust, sep="\t", index=False, header=False)
strchive["motif_length"] = strchive["info"].apply(lambda x: parse_motifs(x))
# get the number of copies of the motif in the reference genome
strchive["num_copies"] = (strchive["end"] - strchive["start"]) / strchive[
    "motif_length"
]
targets_longtr = os.path.join(work_dir, "data/STRchive_LongTR.bed")
strchive[["chrom", "start", "end", "motif_length", "num_copies"]].to_csv(
    targets_longtr,
    sep="\t",
    index=False,
    header=False,
)


if "extra" not in config:
    config["extra"] = ""

# setting up sample dataframes

# VIENNA
samples_vienna = pd.read_table(
    os.path.join(work_dir, "data/all_cram_names_ftp.txt"),
    header=None,
    names=["filename"],
)
samples_vienna_ftp = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/"
)
samples_vienna["sample"] = samples_vienna["filename"].str.replace(
    ".hg38.cram", "", regex=False
)
samples_vienna["hg38_path"] = (
    samples_vienna_ftp + "hg38/" + samples_vienna["filename"].astype(str)
)
samples_vienna["t2t_path"] = (
    samples_vienna_ftp
    + "t2t/"
    + samples_vienna["filename"].str.replace(".hg38.", ".t2t.", regex=False).astype(str)
)
samples_vienna["source"] = "Noyvert/Schloissnig"

# MILLER
samples_miller = pd.read_table(
    os.path.join(work_dir, "data/list-miller-20240606.txt"),
    header=None,
    names=["filename"],
)

samples_miller_s3 = "https://s3.amazonaws.com/1000g-ont/"

# This line is awful.
samples_miller["sample"] = (
    samples_miller["filename"]
    .apply(lambda x: os.path.basename(x))
    .str.split("-")
    .str[0]
    .str.split("_")
    .str[0]
    .str.replace("GM", "NA", regex=False)
)

samples_miller["hg38_path"] = samples_miller_s3 + samples_miller["filename"].astype(str)
samples_miller["t2t_path"] = "not today"
samples_miller["source"] = "Gustafson"


# SAMPLEINFO
sample_info = pd.read_table(
    os.path.join(work_dir, "data/igsr_samples.tsv"),
    usecols=["Sample name", "Superpopulation code", "Sex"],
).set_index("Sample name")

# Samples with missing Sex or Superpopulation code are dropped
# duplicates are dropped, keeping the last one as samples from Miller are generally better covered
samples = (
    pd.concat([samples_vienna, samples_miller])
    .drop_duplicates(subset="sample", keep="last")
    .set_index("sample")
    .join(sample_info)
    .dropna()
)
print("\nSources and samples in study:")
print(samples["source"].value_counts().to_string() + "\n")

samples.to_csv(
    os.path.join(work_dir, "data/pathSTR_samples.tsv"),
    sep="\t",
    index=True,
)


def get_ref(wildcards):
    return {
        "hg38": "/home/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        "t2t": "/home/wdecoster/database/1KG_ONT_VIENNA_t2t.fa",
    }[wildcards.build]


def get_path(wildcards):
    return samples.loc[wildcards.sample, f"{wildcards.build}_path"]


def get_haploid_chroms(wildcards):
    if samples.loc[wildcards.sample, "Sex"] == "female":
        # all chromosomes are diploid, argument should not be used
        return ""
    else:
        return "--haploid chrX,chrY"


def get_sex(wildcards):
    return samples.loc[wildcards.sample, "Sex"]


rule all:
    input:
        strdust=expand(
            os.path.join(work_dir, "pathSTR_STRdust/{build}/{sample}.vcf.gz"),
            sample=samples.index,
            build=["hg38"],
        ),
        longtr=expand(
            os.path.join(work_dir, "pathSTR_longTR/{build}/{sample}.vcf.gz"),
            sample=samples.index,
            build=["hg38"],
        ),
        length_vs_yield=os.path.join(work_dir, "plots/yield_vs_length.html"),
        good_samples=expand(
            os.path.join(
                work_dir, "pathSTR_{genotyper}_good_samples/hg38/good_samples.txt"
            ),
            genotyper=["STRdust", "longTR"],
        ),
        good_samples_zip=expand(
            os.path.join(work_dir, "pathSTR_{genotyper}_good_samples.zip"),
            genotyper=["STRdust", "longTR"],
        ),
        sex_check=os.path.join(work_dir, "plots/sex_check.html"),


rule random_repeats:
    input:
        strdust_random=expand(
            os.path.join(
                work_dir, "pathSTR_STRdust_random_repeats/hg38/{sample}.vcf.gz"
            ),
            sample=samples.index,
        ),


rule extra_repeats:
    # using --config to specify the extra repeats
    input:
        strdust_extra=expand(
            os.path.join(
                work_dir, "pathSTR_STRdust_extra_repeats/hg38/{sample}.vcf.gz"
            ),
            sample=samples.index,
        ),


rule strdust_unphased:
    output:
        os.path.join(work_dir, "pathSTR_STRdust/{build}/{sample}.vcf.gz"),
    log:
        os.path.join(work_dir, "logs/pathSTR_STRdust/{build}/{sample}.log"),
    params:
        cram=get_path,
        ref=get_ref,
        targets=targets_strdust,
        binary=strdust,
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


rule longTR:
    # this rule uses a hacked version of LongTR, in which the check for existance of the input cram is removed (as the file is remote)
    output:
        os.path.join(work_dir, "pathSTR_longTR/{build}/{sample}.vcf.gz"),
    log:
        os.path.join(work_dir, "logs/pathSTR_longTR/{build}/{sample}.log"),
    params:
        cram=get_path,
        ref=get_ref,
        sample="{wildcards.sample}",
        sex=get_sex,
        targets=targets_longtr,
        binary=longTR,
    shell:
        """{params.binary} \
        --bams {params.cram} \
        --fasta {params.ref} \
        --regions {params.targets} \
        --bam-samps {params.sample} \
        --bam-libs lib1 \
        --tr-vcf {output} \
        --min-mean-qual -1 2> {log}
        """


rule strdust_unphased_random_repeats:
    output:
        os.path.join(work_dir, "pathSTR_STRdust_random_repeats/{build}/{sample}.vcf.gz"),
    log:
        os.path.join(
            work_dir, "logs/pathSTR_STRdust_random_repeats/{build}/{sample}.log"
        ),
    params:
        cram=get_path,
        ref=get_ref,
        targets="/home/wdecoster/pathSTR-1000G/data/random-simple-repeats.bed",
        binary=strdust,
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        """RUST_LOG=debug {params.binary} \
        -R {params.targets} \
        --support 2 \
        --unphased \
        --find-outliers \
        --somatic \
        --threads 1 \
        {params.ref} \
        {params.cram} 2> {log} \
        | bgzip > {output} 2>> {log}
        """


rule strdust_unphased_extra_repeats:
    output:
        os.path.join(work_dir, "pathSTR_STRdust_extra_repeats/{build}/{sample}.vcf.gz"),
    log:
        os.path.join(
            work_dir, "logs/pathSTR_STRdust_extra_repeats/{build}/{sample}.log"
        ),
    params:
        cram=get_path,
        ref=get_ref,
        targets=config["extra"],
        binary=strdust,
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        """RUST_LOG=debug {params.binary} \
        -R {params.targets} \
        --support 2 \
        --unphased \
        --find-outliers \
        --somatic \
        --threads 1 \
        {params.ref} \
        {params.cram} 2> {log} \
        | bgzip > {output} 2>> {log}
        """


rule cramino:
    output:
        os.path.join(work_dir, "cramino/{build}/{sample}.cramino"),
    log:
        os.path.join(work_dir, "logs/cramino/{build}/{sample}.log"),
    params:
        cram=get_path,
        ref=get_ref,
        binary="/home/wdecoster/repositories/cramino/target/release/cramino",
    shell:
        "{params.binary} --karyotype --reference {params.ref} {params.cram} > {output} 2> {log}"


rule cramino_gather:
    input:
        expand(
            os.path.join(work_dir, "cramino/hg38/{sample}.cramino"),
            sample=samples.index,
        ),
    output:
        os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
    log:
        "logs/cramino_gather.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(work_dir, "scripts/cramino_gather.py"),
    shell:
        "/home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -i {input} -o {output} 2> {log}"


rule plot_length_vs_yield:
    input:
        os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
    output:
        os.path.join(work_dir, "plots/yield_vs_length.html"),
    log:
        "logs/yield_vs_length.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(work_dir, "scripts/yield_vs_length.py"),
        sample_info=os.path.join(work_dir, "data/pathSTR_samples.tsv"),
    shell:
        "/home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -i {input} -s {params.sample_info} -o {output} 2> {log}"


rule copy_good_samples:
    input:
        overview=os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
        vcfs=expand(
            os.path.join(work_dir, "pathSTR_{{genotyper}}/hg38/{sample}.vcf.gz"),
            sample=samples.index,
        ),
    output:
        os.path.join(work_dir, "pathSTR_{genotyper}_good_samples/hg38/good_samples.txt"),
    log:
        "logs/copy_good_samples_{genotyper}.log",
    conda:
        "/home/wdecoster/pathSTR-1000G/envs/pandas_plotly.yml"
    params:
        script=os.path.join(work_dir, "scripts/copy_good_samples.py"),
    shell:
        "work_dir=$(dirname {output}) ; /home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} -c {input.overview} -v {input.vcfs} -o $work_dir 2> {log}"


rule zip_good_samples:
    input:
        os.path.join(work_dir, "pathSTR_{genotyper}_good_samples/hg38/good_samples.txt"),  #doesn't actually use that file, just needs to know it has been created
    output:
        os.path.join(work_dir, "pathSTR_{genotyper}_good_samples.zip"),
    log:
        "logs/zip_good_samples_{genotyper}.log",
    shell:
        "work_dir=$(dirname {input}) && zip -r -j {output} $work_dir &> {log}"


rule plot_sex_check:
    input:
        os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
    output:
        os.path.join(work_dir, "plots/sex_check.html"),
    log:
        "logs/sex_check.log",
    params:
        script=os.path.join(os.path.dirname(workflow.basedir), "scripts/sex_check.py"),
        sampleinfo=os.path.join(work_dir, "data/pathSTR_samples.tsv"),
        minyield=32,
    shell:
        """/home/wdecoster/miniconda3/envs/pandas_plotly/bin/python {params.script} \
        --cramino {input} \
        --sampleinfo {params.sampleinfo} \
        --minyield {params.minyield} > {output} 2> {log}"""
