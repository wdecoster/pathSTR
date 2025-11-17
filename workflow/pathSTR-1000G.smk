import os
import pandas as pd


def parse_motifs(info):
    """Return the motif of the repeat"""
    return [i.split('=')[1] for i in info.split(";") if i.startswith("MOTIFS=")][0] 

# setting up paths
work_dir = "/home/AD/wdecoster/pathSTR/"
strdust = "/home/AD/wdecoster/repositories/STRdust/target/release/STRdust"

## setting up loci to target from STRchive
strchive = pd.read_csv(
    "https://github.com/dashnowlab/STRchive/raw/main/data/catalogs/STRchive-disease-loci.hg38.TRGT.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "info"],
)
targets_strdust = os.path.join(work_dir, "data/STRchive_STRdust.bed")
strchive.to_csv(targets_strdust, sep="\t", index=False, header=False)
strchive["motif"] = strchive["info"].apply(lambda x: parse_motifs(x))

targets_longtr = os.path.join(work_dir, "data/STRchive_LongTR.bed")
strchive[["chrom", "start", "end", "motif"]].to_csv(
    targets_longtr,
    sep="\t",
    index=False,
    header=False,
)
strchive_t2t = pd.read_csv(
    "https://github.com/dashnowlab/STRchive/raw/main/data/catalogs/STRchive-disease-loci.T2T-chm13.TRGT.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "info"],
)
targets_strdust_t2t = os.path.join(work_dir, "data/STRchive_STRdust_t2t.bed")
strchive_t2t.to_csv(targets_strdust_t2t, sep="\t", index=False, header=False)
strchive_t2t["motif"] = strchive_t2t["info"].apply(lambda x: parse_motifs(x))
targets_longtr_t2t = os.path.join(work_dir, "data/STRchive_LongTR_t2t.bed")
strchive_t2t[["chrom", "start", "end", "motif"]].to_csv(
    targets_longtr_t2t,
    sep="\t",
    index=False,
    header=False,
)


if "extra" not in config:
    config["extra"] = ""

# setting up sample dataframes

def vienna():
    df = pd.read_table(
        os.path.join(work_dir, "data/all_cram_names_ftp.txt"),
        header=None,
        names=["filename"],
    )
    vienna_ftp = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/"
)
    df["sample"] = df["filename"].str.replace(".hg38.cram", "", regex=False)
    
    df["hg38_path"] = (vienna_ftp + "hg38/" + df["filename"].astype(str))
    df["t2t_path"] = (vienna_ftp + "t2t/" + df["filename"].str.replace(".hg38.", ".t2t.", regex=False).astype(str))
    df["source"] = "Noyvert/Schloissnig"
    return df.drop(columns="filename")

def gustafson(build):
    df = pd.read_table(
        os.path.join(
            work_dir, f"data/list-miller-2025103_{build}.txt"
        ),
        header=None,
        names=["filename"],
    )
    df["sample"] = (
        df["filename"]
        .apply(lambda x: os.path.basename(x))
        .str.split("-")
        .str[0]
        .str.split("_")
        .str[0]
        .str.replace("GM", "NA", regex=False)
    )
    buildname = {"hg38": "hg38", "chm13": "t2t"}[build]
    df[f"{buildname}_path"] = (
        "https://s3.amazonaws.com/1000g-ont/" + df["filename"].astype(str)
    )
    df["source"] = "Gustafson"
    return df.drop(columns="filename")

# GUSTAFSON
# note that this function contains hardcoded paths to a list generated using aws s3 ls, which will have to be updated when updating the database
samples_vienna = vienna()
samples_gustafson_hg38 = gustafson("hg38")
samples_gustafson_t2t = gustafson("chm13")
# both tables from gustafson have sample, hg38_path/t2t_path and source columns, so can be merged directly
samples_gustafson = pd.merge(
    samples_gustafson_hg38,
    samples_gustafson_t2t,
    on=["sample", "source"],
    how="outer",
)

# SAMPLEINFO
sample_info = pd.read_table(
    os.path.join(work_dir, "data/igsr_samples.tsv"),
    usecols=["Sample name", "Superpopulation code", "Sex"],
).set_index("Sample name")

# Samples with missing Sex or Superpopulation code are dropped
# duplicates are dropped, keeping the last one as samples from Miller are generally better covered
samples = (
    pd.concat([samples_vienna, samples_gustafson])
    .drop_duplicates(subset="sample", keep="last")
    .set_index("sample")
    .join(sample_info)
    .dropna(subset=["Sex", "Superpopulation code"])
    
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
        "hg38": "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        "t2t": "/home/AD/wdecoster/database/1KG_ONT_VIENNA_t2t.fa",
    }[wildcards.build]


def get_targets_strdust(wildcards):
    return {
        "hg38": targets_strdust,
        "t2t": targets_strdust_t2t,
    }[wildcards.build]


def get_targets_longtr(wildcards):
    return {
        "hg38": targets_longtr,
        "t2t": targets_longtr_t2t,
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


def get_vcfs(wildcards):
    return [
        os.path.join(
            work_dir, f"pathSTR_{wildcards.genotyper}/{wildcards.build}/{s}.vcf.gz"
        )
        for s in samples.dropna(subset=[f"{wildcards.build}_path"]).index.to_list()
    ]


rule all:
    input:
        strdust=expand(
            os.path.join(work_dir, "pathSTR_STRdust/{build}/{sample}.vcf.gz"),
            sample=samples.dropna(subset=["hg38_path"]).index,
            build=["hg38"],
        ),
        longtr=expand(
            os.path.join(work_dir, "pathSTR_LongTR/{build}/{sample}.vcf.gz"),
            sample=samples.dropna(subset=["hg38_path"]).index,
            build=["hg38"],
        ),
        length_vs_yield=os.path.join(work_dir, "plots/yield_vs_length.html"),
        good_samples=expand(
            os.path.join(
                work_dir, "pathSTR_{genotyper}_good_samples/{build}/good_samples.txt"
            ),
            genotyper=["STRdust", "LongTR"],
            build=["hg38", "t2t"],
        ),
        good_samples_zip=expand(
            os.path.join(work_dir, "pathSTR_{genotyper}_{build}_good_samples.zip"),
            genotyper=["STRdust", "LongTR"],
            build=["hg38", "t2t"],
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
        targets=get_targets_strdust,
        binary=strdust,
        haploid_chroms=get_haploid_chroms,
    conda:
        "/home/AD/wdecoster/anaconda3/envs/samtools"
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
        os.path.join(work_dir, "pathSTR_LongTR/{build}/{sample}.vcf.gz"),
    log:
        os.path.join(work_dir, "logs/pathSTR_LongTR/{build}/{sample}.log"),
    params:
        cram=get_path,
        ref=get_ref,
        sample=lambda wildcards: {wildcards.sample},
        sex=get_sex,
        targets=get_targets_longtr,
        binary="~/anaconda3/envs/longtr/bin/LongTR"
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
        "/home/AD/wdecoster/anaconda3/envs/samtools"
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
        "/home/AD/wdecoster/anaconda3/envs/samtools"
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
        binary="/home/AD/wdecoster/repositories/cramino/target/release/cramino",
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
        os.path.join(work_dir, "envs/pandas_plotly.yml")
    params:
        script=os.path.join(work_dir, "scripts/cramino_gather.py"),
    shell:
        "python {params.script} -i {input} -o {output} 2> {log}"


rule plot_length_vs_yield:
    input:
        os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
    output:
        os.path.join(work_dir, "plots/yield_vs_length.html"),
    log:
        "logs/yield_vs_length.log",
    conda:
        os.path.join(work_dir, "envs/pandas_plotly.yml")
    params:
        script=os.path.join(work_dir, "scripts/yield_vs_length.py"),
        sample_info=os.path.join(work_dir, "data/pathSTR_samples.tsv"),
    shell:
        "python {params.script} -i {input} -s {params.sample_info} -o {output} 2> {log}"


# this rule uses, regardless of the build, the hg38 good samples
# this is under the assumption that the good samples are the same for both builds
# and that for every sample with a t2t path, there is also a hg38 path
# the other way around is not necessarily true
# this may eventually be a problem, but for now it avoids to rerun cramino on the same samples for another build
rule copy_good_samples:
    input:
        vcfs=get_vcfs,
        overview=os.path.join(work_dir, "cramino/hg38/cramino_all.tsv"),
    output:
        os.path.join(
            work_dir, "pathSTR_{genotyper}_good_samples/{build}/good_samples.txt"
        ),
    log:
        "logs/copy_good_samples_{genotyper}_{build}.log",
    conda:
        os.path.join(work_dir, "envs/pandas_plotly.yml")
    params:
        script=os.path.join(work_dir, "scripts/copy_good_samples.py"),
    shell:
        "work_dir=$(dirname {output}) ; python {params.script} -c {input.overview} -v {input.vcfs} -o $work_dir 2> {log}"


rule zip_good_samples:
    input:
        os.path.join(
            work_dir, "pathSTR_{genotyper}_good_samples/{build}/good_samples.txt"
        ),
        #doesn't actually use that file, just needs to know it has been created
    output:
        os.path.join(work_dir, "pathSTR_{genotyper}_{build}_good_samples.zip"),
    log:
        "logs/zip_good_samples_{genotyper}_{build}.log",
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
    conda:
        os.path.join(work_dir, "envs/pandas_plotly.yml")
    shell:
        """python {params.script} \
        --cramino {input} \
        --sampleinfo {params.sampleinfo} \
        --minyield {params.minyield} > {output} 2> {log}"""
