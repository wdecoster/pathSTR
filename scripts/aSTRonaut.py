# script to create a repeat sequence plot similar to the pathSTR <sequence> composition vizualization, but stand-alone
# taking VCF files as input

from argparse import ArgumentParser
import plotly.express as px
import pandas as pd
from collections import Counter
from cyvcf2 import VCF
from itertools import chain
import os

def main():
    args = get_args()
    df = parse_input(args.vcf, args.repeat)
    with open(args.out, "w") as out:
        for repeat in df["coords"].unique():
            repeat_df = df[df["coords"] == repeat]
            kmer_df = parse_kmers(repeat_df, args.kmer)
            plot_sequence(repeat_df, kmer_df, repeat).write_html(out)        

def plot_sequence(repeat_df, kmer_df, repeat):
    """
    Plot the sequence composition of a repeat
    assign a color to each kmer
    using only the 10 most frequent kmers by limiting kmer_df to the first 10 columns
    :param repeat_df: dataframe with the repeat sequences
    :param kmer_df: dataframe with the kmer counts
    :param repeat: coordinates of the repeat
    """
    colors = [
        "rgb(31, 119, 180)",
        "rgb(255, 127, 14)",
        "rgb(44, 160, 44)",
        "rgb(214, 39, 40)",
        "rgb(148, 103, 189)",
        "rgb(140, 86, 75)",
        "rgb(227, 119, 194)",
        "rgb(127, 127, 127)",
        "rgb(188, 189, 34)",
        "rgb(23, 190, 207)",
    ]
    kmer_dict = {k: c for k, c in zip(kmer_df.columns[:10], colors)}
    inverse_dict = {v: k for k, v in kmer_dict.items()}
    inverse_dict["rgb(128, 128, 128)"] = "other"
    repeat_df = repeat_df.dropna(subset=["sequence"]).sort_values(by="sequence", key = lambda x: x.str.len())
    # draw a scatter plot showing for each sample and allele the order of the kmers in the sequence
    # replace in the sequence from frequent to less frequent the kmers with their color
    repeat_df["seq_colored"] = repeat_df["sequence"]
    for k, c in kmer_dict.items():
        repeat_df["seq_colored"] = repeat_df["seq_colored"].str.replace(
            k, f"{c};" * len(k)
        )
    # replace the remaining nucleotides with grey
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.replace(
        r"[ACTG]", "rgb(128, 128, 128);", regex=True
    )
    # remove the last semicolon
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.rstrip(";")
    # split the seq_colored into a list
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.split(";")
    # add a range column to the dataframe, to enumerate the nucleotides
    repeat_df["range"] = repeat_df["seq_colored"].apply(
        lambda x: list(range(len(x)))
    )
    repeat_df["identifier"] = repeat_df["sample"] + "_" + repeat_df["allele"]
    # explode the seq_colored and range columns for plotting
    repeat_colors = repeat_df[
        ["identifier", "sequence", "seq_colored", "range"]
    ].explode(["seq_colored", "range"])
    # convert colors back to kmers
    repeat_colors["kmer"] = repeat_colors["seq_colored"].apply(
        lambda x: inverse_dict[x]
    )

    fig = px.scatter(
        repeat_colors,
        x="range",
        y="identifier",
        color="kmer",
        color_discrete_sequence=repeat_colors["seq_colored"].unique(),
        hover_data=["kmer"],
        labels={
            "range": "Nucleotide position in repeat",
            "identifier": "Sample/allele",
        },
        title=f"Repeat {repeat}",
        category_orders={"identifier": repeat_df["identifier"][::-1]},
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_yaxes(tickfont_size=8)
    return fig

def parse_kmers(repeat_df, motif_length):
    kmers_extracted = []
    for sample, allele, seq in repeat_df[["sample", "allele", "sequence"]].itertuples(
        index=False, name=None
    ):
        if seq:
            kmers = count_kmers(seq, k=motif_length)
            if kmers:
                kmers.update(
                    {
                        "identifier": f"{sample}_{allele}",
                        "length": len(seq) / motif_length,
                    }
                )
                kmers_extracted.append(kmers),
    return (
        pd.DataFrame(kmers_extracted)
        .set_index("identifier")
        .astype(float)
        .fillna(0.0)
        .round(2)
    )


def count_kmers(seq, k):
    """
    Count the kmers in a sequence
    :param seq: sequence to count kmers in
    :param k: kmer length
    """
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmers[seq[i : i + k]] += 1
    return prune_counts(kmers)


def get_rotations(kmer):
    """
    Rotate a kmer to get all equivalent representations
    Return the lexicographical first separately. Also return all rotations.
    """
    e = len(kmer)
    rotations = [kmer[i:e] + kmer[:i] for i in range(e)]
    return sorted(rotations)[0], rotations


def prune_counts(kmers):
    """
    For all rotations of a kmer, keep only the lexicographical first
    Return the number as a fraction of the total kmers
    """
    pruned = dict()
    for key in kmers:
        first, all_rotations = get_rotations(key)
        if first in pruned.keys():
            continue
        else:
            pruned[first] = sum([kmers[r] for r in all_rotations])
    total_kmers = sum(pruned.values())
    return {k: v / total_kmers for k, v in pruned.items()}

def parse_input(vcf_list, repeat):
    """
    Parse the VCF files and return a dataframe
    :param vcf_list: list of VCF files to parse
    :param repeat: coordinates of the repeat to extract, or None if all have to be extracted
    """
    # read in the VCFs
    calls = [parse_vcf(vcf, repeat) for vcf in vcf_list]
    # make a dataframe and join with the sample info
    df = pd.DataFrame(
            flatten(calls),
            columns=[
                "coords",
                "sample",
                "allele",
                "sequence",
            ],
        ).set_index("sample", drop=False)
    return df


def parse_vcf(vcf, repeat=None):
    """
    Parse a VCF file and return a list of sequences
    :param vcf: path to the VCF file
    :param repeat: coordinates of the repeat to extract, or None if all have to be extracted
    """
    calls = []
    name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        coords = f"{v.CHROM}:{v.POS}"
        if repeat and coords != repeat:
            continue
        sequences = parse_alts(v.ALT, v.genotypes[0])
        calls.append((coords, name, "Allele1", sequences[0]))
        calls.append((coords, name, "Allele2", sequences[1]))
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



def get_args():
    parser = ArgumentParser(description="Create a repeat sequence plot similar to the pathSTR <sequence> composition vizualization, but stand-alone")
    parser.add_argument("vcf", help="VCF files to analyze", nargs="+")
    parser.add_argument("-k", "--kmer", help="Kmer length to use for plot", default=3, type=int)
    parser.add_argument("--repeat", help="Chromosome and POS of repeat to plot e.g. chr1:12345", default=None)
    parser.add_argument("-o", "--out", help="Output file name", default="reptor.html")
    return parser.parse_args()

if __name__ == "__main__":
    main()