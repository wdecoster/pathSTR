import pandas as pd
from itertools import chain
import os
from cyvcf2 import VCF
import gzip
import base64
import sys


def parse_input(vcf_list, sample_info, repeats):

    # read in the VCFs
    lengths = [get_lengths_from_vcf(vcf, repeats) for vcf in vcf_list]
    # make a dataframe and join with the sample info
    df = (
        pd.DataFrame(
            flatten(lengths),
            columns=[
                "chrom",
                "gene",
                "sample",
                "allele",
                "length",
                "ref_diff",
                "sequence",
                "support",
            ],
        )
        .set_index("sample", drop=False)
        .join(
            pd.read_csv(
                sample_info,
                sep="\t",
                usecols=["Sample name", "Sex", "Superpopulation code"],
            )
            .set_index("Sample name")
            .rename(columns={"Superpopulation code": "Superpopulation"})
        )
    ).assign(Group="1000 Genomes")
    # for every repeat in the dataframe, divide the length by the motif length
    df["length"] = df.apply(
        lambda x: x["length"] / repeats.motif_length(x["gene"]), axis=1
    )
    df["ref_diff"] = df.apply(
        lambda x: x["ref_diff"] / repeats.motif_length(x["gene"]), axis=1
    )
    # Remove duplicate hits for males on chrX
    df = df[
        ~(
            (df["chrom"] == "chrX")
            & (df["Sex"] == "male")
            & df.duplicated(subset=["sample", "gene", "length"], keep="first")
        )
    ]
    return df


def get_lengths_from_vcf(vcf, repeats):
    calls = []
    name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        gene = repeats.gene(f"{v.CHROM}:{str(v.POS)}-{str(v.end)}")
        if gene is None:
            print(f"Skipping {v.CHROM}:{str(v.POS)}-{str(v.end)} - not in bed file.")
            continue
        full_lengths = v.INFO.get("FRB")
        ref_diff = v.INFO.get("RB")
        support = v.format("SUP")[0]
        sequences = parse_alts(v.ALT, v.genotypes[0])
        calls.append(
            (
                v.CHROM,
                gene,
                name,
                "Allele1",
                full_lengths[0],
                ref_diff[0],
                sequences[0],
                support[0],
            )
        )
        calls.append(
            (
                v.CHROM,
                gene,
                name,
                "Allele2",
                full_lengths[1],
                ref_diff[1],
                sequences[1],
                support[1],
            )
        )
    return calls


def parse_alts(alts, genotype):
    sequences = []
    for phase in [0, 1]:
        if genotype[phase] in [0, -1]:
            sequences.append(None)
        elif genotype[phase] == 1:
            sequences.append(alts[0])
        elif genotype[phase] == 2:
            sequences.append(alts[1])
        else:
            raise ValueError("Unexpected genotype")
    return sequences


def flatten(it):
    return chain.from_iterable(it)


def get_lengths_from_uploaded_vcf(contents, filename, repeats):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    # write the decoded file to a temporary file
    # make sure the file ends with .gz and is gzipped
    # and read it back in using cyvcf2
    if filename.endswith(".gz"):
        tempfile = os.path.join("/tmp", os.path.basename(filename))
        with open(tempfile, "wb") as f:
            f.write(decoded)
    else:
        tempfile = os.path.join("/tmp", os.path.basename(filename) + ".gz")
        with gzip.open(tempfile, "wb") as f:
            f.write(decoded)
    try:
        calls = get_lengths_from_vcf(tempfile, repeats)
    except OSError:
        # this happens when the file is not a VCF, or malformed
        os.remove(tempfile)
        return None
    df = pd.DataFrame(
        calls,
        columns=[
            "chrom",
            "gene",
            "sample",
            "allele",
            "length",
            "ref_diff",
            "sequence",
            "support",
        ],
    )
    # for every repeat in the dataframe, divide the length and ref_diff by the motif length
    df["length"] = df.apply(
        lambda x: round(x["length"] / repeats.motif_length(x["gene"]), 2),
        axis=1,
    )
    df["ref_diff"] = df.apply(
        lambda x: round(x["ref_diff"] / repeats.motif_length(x["gene"]), 2),
        axis=1,
    )
    df["Group"] = "Uploaded"
    df["Superpopulation"] = "Uploaded"
    df["Sex"] = "Uploaded"
    os.remove(tempfile)
    return df


def stats(df):
    """Calculate summary statistics mean and standard deviation for the data"""
    return df.groupby("gene")["length"].agg(["mean", "std"]).round(1)


def get_composition(df, gene, repeats):
    """Calculate the composition of the dataset"""
    motifs = repeats.motifs(gene)
    df = df[df["gene"] == gene].drop(columns=["chrom", "length", "ref_diff"])
