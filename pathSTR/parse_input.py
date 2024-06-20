import pandas as pd
from itertools import chain
import os
import io
import gzip
import base64
from cyvcf2 import VCF
import numpy as np
from .rle import rle
import logging


def parse_input(vcf_fofn, sample_info, repeats):
    # read the file of filenames
    fofn = pd.read_csv(
        vcf_fofn, header=None, names=["vcf", "build", "caller"], sep="\t"
    )
    # read in the VCFs
    calls = [
        parse_vcf(vcf, build, caller, repeats)
        for (vcf, build, caller) in fofn.itertuples(index=False, name=None)
    ]
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
                "dataset",
            ],
        )
        .set_index("sample", drop=False)
        .join(
            pd.read_csv(
                sample_info,
                sep="\t",
                usecols=[
                    "sample",
                    "Sex",
                    "Superpopulation code",
                    "source",
                    "hg38_path",
                ],
            )
            .set_index("sample")
            .rename(columns={"Superpopulation code": "Superpopulation"})
        )
    ).assign(Group="1000 Genomes")

    # for every repeat in the dataframe, divide the length by the motif length
    # playing a bit of a tricky game here, as the motif length for hg38 is used
    # this is under the very reasonable assumption that the motif length is the same for all genome builds
    # it is, at this stage, rather hard to run this function for multiple builds
    # but maybe one day this comes back to bite us
    df["length"] = df.apply(
        lambda x: x["length"] / repeats.motif_length(x["gene"], build="hg38"), axis=1
    ).round(2)
    df["ref_diff"] = df.apply(
        lambda x: x["ref_diff"] / repeats.motif_length(x["gene"], build="hg38"), axis=1
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


def parse_vcf(vcf, build, caller, repeats, name=None):
    """
    Parse a VCF file and return a list of calls
    The chromosome and position are used to get the gene from the repeats object
    If necessary, the chromosome is prefixed with "chr"
    Positions are matched exactly, so the bed file for genotyping should be the same as the one used for the repeats object here

    This function expects that the VCF contains only one sample
    """
    caller_ = caller.lower()
    calls = []
    if name is None:
        name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        chrom = v.CHROM if v.CHROM.startswith("chr") else "chr" + v.CHROM
        if caller_ == "strdust":
            gene = repeats.coords_to_gene(
                f"{chrom}:{str(v.POS)}-{str(v.end)}", build=build
            )
        elif caller_ == "longtr":
            # LongTR may adjust the POS field to include a SNV https://github.com/gymrek-lab/LongTR/issues/8
            start = v.INFO.get("START")
            gene = repeats.coords_to_gene(
                f"{chrom}:{str(start)}-{str(v.end)}", build=build
            )
        else:
            raise ValueError("Unexpected caller")
        if gene is None:
            logging.warning(
                f"Skipping {chrom}:{str(v.POS)}-{str(v.end)} from {caller}:{build} in {vcf} - interval not found in bed file."
            )
            continue
        if caller_ == "strdust":
            full_lengths = v.format("FRB")[0]
            ref_diff = v.format("RB")[0]
            support = v.format("SUP")[0]
        elif caller_ == "longtr":
            # GB is presented as 'x|y' with the difference with the reference per allele
            ref_diff = [int(i) for i in v.format("GB")[0].split("|")]
            # LongTR may adjust the POS field to include a SNV, and does not have a field for the full repeat length
            full_lengths = (
                ref_diff[0] + v.end - start + 1,
                ref_diff[1] + v.end - start + 1,
            )
            # LongTR does not have a field for the support, so I use the depth, which is the total number of reads
            # this is however not the support per allele, so I duplicate the total depth for both alleles
            support = (f"total:{v.format('DP')[0][0]}", f"total:{v.format('DP')[0][0]}")
        else:
            raise ValueError("Unexpected caller")
        sequences = parse_alts(v.ALT, v.genotypes[0])
        # missing genotypes end up as -2147483648 for FRB and RB.
        # that is odd and problematic, so I replace them with np.nan
        # I don't know if that is only the case for STRdust, or if it is a general issue
        calls.append(
            (
                chrom,
                gene,
                name,
                "Allele1",
                full_lengths[0] if full_lengths[0] > -2147483648 else np.nan,
                ref_diff[0] if ref_diff[0] > -2147483648 else np.nan,
                sequences[0],
                support[0],
                f"{caller}_{build}",
            )
        )
        calls.append(
            (
                v.CHROM,
                gene,
                name,
                "Allele2",
                full_lengths[1] if full_lengths[1] > -2147483648 else np.nan,
                ref_diff[1] if ref_diff[1] > -2147483648 else np.nan,
                sequences[1],
                support[1],
                f"{caller}_{build}",
            )
        )
    return calls


def parse_alts(alts, genotype):
    sequences = []
    for phase in [0, 1]:
        if genotype[phase] in [0, -1]:
            sequences.append(None)
        elif genotype[phase] == 1:
            sequences.append(alts[0] if not alts[0] == "<DEL>" else None)
        elif genotype[phase] == 2:
            sequences.append(alts[1] if not alts[1] == "<DEL>" else None)
        else:
            raise ValueError("Unexpected genotype")
    return sequences


def flatten(it):
    return chain.from_iterable(it)


def parse_uploaded_vcf(contents, uploaded_filename, repeats, build, caller):
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
        logging.info(f"Decoding {uploaded_filename}")
        r, w = os.pipe()
        name = uploaded_filename.replace(".vcf", "").replace(".gz", "")
        if uploaded_filename.endswith(".gz"):
            with os.fdopen(w, "wb") as f:
                f.write(decoded)
        else:
            with gzip.open(w, "wb") as f:
                f.write(decoded)
        r_fd = os.fdopen(r, "rb")
        calls = parse_vcf(r_fd, build, caller, repeats, name=name)
        logging.info(f"Parsed {uploaded_filename}")
    except OSError as e:
        # this happens when the file is not a VCF, or malformed
        logging.warning(e)
        return None
    except Exception as e:
        logging.error(e)
        return None
    else:
        if len(calls) == 0:
            logging.warning(f"No valid calls found in {uploaded_filename}")
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
                "dataset",
            ],
        )
        # for every repeat in the dataframe, divide the length and ref_diff by the motif length
        df["length"] = df.apply(
            lambda x: round(
                x["length"] / repeats.motif_length(x["gene"], build=build), 2
            ),
            axis=1,
        )
        df["ref_diff"] = df.apply(
            lambda x: round(
                x["ref_diff"] / repeats.motif_length(x["gene"], build=build), 2
            ),
            axis=1,
        )
        df["Group"] = "Uploaded"
        df["Superpopulation"] = "Uploaded"
        df["Sex"] = "Uploaded"
        return df


def stats(df):
    """Calculate summary statistics mean and standard deviation for the data"""
    return df.groupby(["dataset", "gene"])["length"].agg(["mean", "std"]).round(1)


def create_details_table(df, repeats):
    """
    Reformat the df to an easier format for querying the details per individual
    """
    # playing a bit of a tricky game here, as the motif length for hg38 is used
    # this is under the very reasonable assumption that the motif length is the same for all genome builds
    # it is, at this stage, rather hard to run this function for multiple builds
    # but maybe one day this comes back to bite us
    df["sequence_rle"] = df.apply(
        lambda x: rle(x["sequence"], repeats.motif_length(x["gene"], build="hg38")),
        axis=1,
    )
    logging.info("Creating run-length encoded sequences")
    df = (
        df.drop(columns=["Group", "chrom", "hg38_path", "sequence"])
        .round(1)
        .pivot(
            index=["dataset", "gene", "sample"],
            columns="allele",
            values=[
                "length",
                "ref_diff",
                "sequence_rle",
                "support",
                "Sex",
                "Superpopulation",
                "source",
            ],
        )
        .reset_index()
    )

    # fill in missing columns with the information from the other allele and assign to a new column
    df["sex"] = df[("Sex", "Allele1")].fillna(df[("Sex", "Allele2")])
    df["superpopulation"] = df[("Superpopulation", "Allele1")].fillna(
        df[("Superpopulation", "Allele2")]
    )
    df["data_source"] = df[("source", "Allele2")].fillna(df[("source", "Allele1")])

    return df.drop(
        columns=[
            ("Sex", "Allele1"),
            ("Sex", "Allele2"),
            ("Superpopulation", "Allele1"),
            ("Superpopulation", "Allele2"),
            ("source", "Allele1"),
            ("source", "Allele2"),
        ]
    ).set_index("sample")
