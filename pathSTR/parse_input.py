import pandas as pd
from itertools import chain
import os
import io
import gzip
import base64
from cyvcf2 import VCF


def parse_input(vcf_list, sample_info, repeats):

    # read in the VCFs
    calls = [parse_vcf(vcf, repeats) for vcf in vcf_list]
    # make a dataframe and join with the sample info
    df = (
        pd.DataFrame(
            flatten(calls),
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
    ).round(2)
    df["ref_diff"] = df.apply(
        lambda x: x["ref_diff"] / repeats.motif_length(x["gene"]), axis=1
    ).round(2)
    # Remove duplicate hits for males on chrX
    df = df[
        ~(
            (df["chrom"] == "chrX")
            & (df["Sex"] == "male")
            & df.duplicated(subset=["sample", "gene", "length"], keep="first")
        )
    ]
    return df


def parse_vcf(vcf, repeats, name=None):
    """
    Parse a VCF file and return a list of calls
    The chromosome and position are used to get the gene from the repeats object
    If necessary, the chromosome is prefixed with "chr"
    Positions are matched exactly, so the bed file for genotyping should be the same as the one used for the repeats object here

    This function expects that the VCF contains only one sample
    """
    calls = []
    if name is None:
        name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        chrom = v.CHROM if v.CHROM.startswith("chr") else "chr" + v.CHROM
        gene = repeats.gene(f"{chrom}:{str(v.POS)}-{str(v.end)}")
        if gene is None:
            print(f"Skipping {chrom}:{str(v.POS)}-{str(v.end)} - not in bed file.")
            continue
        full_lengths = v.INFO.get("FRB") if v.INFO.get("FRB") else v.format("FRB")[0]
        ref_diff = v.INFO.get("RB") if v.INFO.get("RB") else v.format("RB")[0]
        support = v.format("SUP")[0]
        sequences = parse_alts(v.ALT, v.genotypes[0])
        calls.append(
            (
                chrom,
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


def parse_uploaded_vcf(contents, uploaded_filename, repeats):
    """
    Parse the uploaded VCF file and return a dataframe
    :param contents: contents of the uploaded file
    :param uploaded_filename: name of the uploaded file
    :param repeats: repeats object

    This function expects that the VCF contains only one sample
    The file is expected to be a VCF or a gzipped VCF, and is decoded from base64, kept in memory, and passed to cyvcf2
    """
    _, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    try:
        r, w = os.pipe()
        name = uploaded_filename.replace(".vcf", "").replace(".gz", "")
        if uploaded_filename.endswith(".gz"):
            with os.fdopen(w, "wb") as f:
                f.write(decoded)
        else:
            with gzip.open(w, "wb") as f:
                f.write(decoded)
        r_fd = os.fdopen(r, "rb")
        calls = parse_vcf(r_fd, repeats, name=name)
    except OSError:
        # this happens when the file is not a VCF, or malformed
        return None
    else:
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
        return df


def stats(df):
    """Calculate summary statistics mean and standard deviation for the data"""
    return df.groupby("gene")["length"].agg(["mean", "std"]).round(1)


# def get_composition(df, gene, repeats):
#     """Calculate the composition of the dataset"""
#     motifs = repeats.motifs(gene)
#     df = df[df["gene"] == gene].drop(columns=["chrom", "length", "ref_diff"])
