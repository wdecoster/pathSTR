# script to create a repeat sequence plot similar to the pathSTR <sequence> composition vizualization, but stand-alone
# taking VCF files as input

from argparse import ArgumentParser
import plotly.express as px
import pandas as pd
from collections import Counter
from cyvcf2 import VCF
from itertools import chain
import os
import sys


def main():
    args = get_args()
    df = parse_input(args)
    if args.sampleinfo:
        sampleinfo = (
            pd.read_table(args.sampleinfo, usecols=["name", "group"])
            .rename(columns={"name": "sample"})
            .set_index("sample")
        )
        # verify that 'case' is a value in the 'group' column
        if "case" not in sampleinfo["group"].unique():
            sys.stderr.write(
                "ERROR: 'case' is not a value in the 'group' column of the sample info file!\n"
            )
            sys.exit(1)
        df = df.join(sampleinfo, how="left")
        # add a 'case' column which is 1 if the sample is a case and 0 otherwise
        df["case"] = df["group"].apply(lambda x: 1 if x == "case" else 0)
    else:
        # if no sample info is provided, all samples are considered controls
        # and no annotation is added to the plot
        df["case"] = 0
    with open(args.out, "w") as out:
        for repeat in df["coords"].unique():
            repeat_df = df[df["coords"] == repeat]
            # Either use the motifs as specified by the user, or select the <number> most frequent kmers
            # with args.motifs the kmer length is ignored
            # if kmers are specified as argument, the kmers are sorted by length
            # so the longest kmers are plotted first, to accomodate for overlapping kmers
            kmers = (
                sorted(args.motifs.split(","), key=lambda x: len(x), reverse=True)
                if args.motifs
                else select_kmers(repeat_df, args.kmer, args.number)
            )
            plot_sequence(repeat_df, kmers, repeat, args).write_html(out)


def plot_sequence(repeat_df, kmers, repeat, args):
    """
    Plot the sequence composition of a repeat
    assign a color to each kmer
    using only the <number> most frequent kmers by limiting kmer_df to the first <number> columns, assuming it is sorted by column sum
    :param repeat_df: dataframe with the repeat sequences
    :param kmers: list of kmers to plot
    :param repeat: coordinates of the repeat
    """
    if len(kmers) > 10:
        colors = [hex_to_rgb(c) for c in px.colors.qualitative.Light24]
        if len(kmers) > len(colors):
            sys.stderr.write("WARNING: Not enough colors defined!\n\n")
    else:
        colors = px.colors.qualitative.Safe
    kmer_dict = {k: c for k, c in zip(kmers, colors)}
    inverse_dict = {v: k for k, v in kmer_dict.items()}
    inverse_dict["rgb(128, 128, 128)"] = "other"

    if args.alphabetic:
        repeat_df["_sample"] = repeat_df["sample"]
        repeat_df = (
            repeat_df.dropna(subset=["sequence"])
            .sort_values(by="_sample", ascending=False)
            .drop(columns="_sample")
        )
    else:
        repeat_df = repeat_df.dropna(subset=["sequence"]).sort_values(
            by="sequence", key=lambda x: x.str.len()
        )
    # draw a scatter plot showing for each sample and allele the order of the kmers in the sequence
    # replace in the sequence the kmers with their color
    repeat_df["seq_colored"] = repeat_df["sequence"]
    for k, c in kmer_dict.items():
        repeat_df["seq_colored"] = repeat_df["seq_colored"].str.replace(
            k, f"{c};" * len(k)
        )
    # replace the remaining nucleotides with black
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.replace(
        r"[ACTG]", "rgb(128, 128, 128);", regex=True
    )
    # remove the last semicolon
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.rstrip(";")
    # split the seq_colored into a list
    repeat_df["seq_colored"] = repeat_df["seq_colored"].str.split(";")
    # add a range column to the dataframe, to enumerate the nucleotides
    repeat_df["range"] = repeat_df["seq_colored"].apply(lambda x: list(range(len(x))))
    if args.hide_allele_label:
        repeat_df["identifier"] = repeat_df["sample"]
    else:
        repeat_df["identifier"] = repeat_df["sample"] + "_" + repeat_df["allele"]
    # explode the seq_colored and range columns for plotting
    repeat_colors = repeat_df[
        ["identifier", "sequence", "seq_colored", "range"]
    ].explode(["seq_colored", "range"])
    # convert colors back to kmers, for the legend and hoverdata
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
            "identifier": "Sample/allele" if not args.somatic else "Reads",
        },
        title=f"Repeat at {repeat}",
        category_orders={"identifier": repeat_df["identifier"][::-1]},
    )

    fig.update_traces(marker=dict(size=args.size, symbol="square"))
    fig.update_xaxes(tickfont_size=20)
    if args.hide_labels:
        fig.update_yaxes(showticklabels=False)
    else:
        fig.update_yaxes(tickfont_size=args.label_size)
    fig.update_layout(legend_traceorder="reversed")

    # add an annotation in front of the case samples, only with --sampleinfo
    for _, row in repeat_df.iterrows():
        if row["case"]:
            fig.add_annotation(
                x=0,
                y=row["identifier"],
                text="<b>>>></b>",
                showarrow=False,
                xshift=-20,
                font=dict(size=10, color="black"),
            )
    if args.publication:
        fig.update_layout(
            plot_bgcolor="rgba(0, 0, 0, 0)",
            paper_bgcolor="rgba(0, 0, 0, 0)",
            font=dict(size=24),
            title=args.title,
        )
        fig.update_xaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
        fig.update_yaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
    fig.update_layout(height=args.height)
    if args.legend_corner == "topright":
        fig.update_layout(
            legend=dict(
                yanchor="top",
                y=0.95,
                xanchor="right",
                x=0.99,
                itemsizing="constant",
            )
        )
    else:
        fig.update_layout(
            legend=dict(
                yanchor="bottom",
                y=0.05,
                xanchor="right",
                x=0.95,
                itemsizing="constant",
            )
        )
    return fig


def hex_to_rgb(hex):
    hex = hex.lstrip("#")
    tup = tuple(int(hex[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgb{tup}"


def select_kmers(repeat_df, motif_length, number=10):
    """
    Function to select the <number> most frequent kmers from the repeat sequences
    :param repeat_df: dataframe with the repeat sequences
    :param motif_length: length of the kmers to count
    :param number: number of kmers to plot
    """
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
    kmer_df = (
        pd.DataFrame(kmers_extracted)
        .set_index("identifier")
        .astype(float)
        .fillna(0.0)
        .round(2)
    )
    # sort the kmer_dfs columns by column sum
    kmer_df = kmer_df[kmer_df.sum().sort_values(ascending=False).index]
    # ignoring the 'length' column, return the <number> most frequent kmers
    return [k for k in kmer_df.columns if k != "length"][:number]


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


def parse_input(args):
    """
    Parse the VCF files and return a dataframe, containing the coordinates, sample name, allele and repeat sequence
    :param args.vcf: list of VCF files to parse
    :param args.repeat: coordinates of the repeat to extract, or None if all have to be extracted
    """
    if args.table:
        df = pd.read_table(args.table)
        if not "name" in df.columns:
            sys.exit("ERROR: No 'name' column found in table")
        if not "sequence" in df.columns:
            sys.exit("ERROR: No 'sequence' column found in table")
        if "group" in df.columns:
            if "case" not in df["group"].unique():
                sys.exit(
                    "ERROR: 'case' is not a value in the 'group' column of the table!\n"
                )
            df["case"] = df["group"].apply(lambda x: 1 if x == "case" else 0)
        else:
            # if no sample info is provided, all samples are considered controls
            # and no annotation is added to the plot
            df["case"] = 0
        df["length"] = df["sequence"].apply(len)
        df = (
            df.loc[df["length"] > args.minlen, ["name", "sequence", "case"]]
            .rename(columns={"name": "sample"})
            .assign(coords="", allele="Allele1")
            .set_index("sample", drop=False)
        )

    else:
        # read in the VCFs
        if args.names:
            calls = [
                parse_vcf(vcf, args, name=name)
                for vcf, name in zip(args.vcf, args.names.split(","))
            ]
        else:
            calls = [parse_vcf(vcf, args) for vcf in args.vcf]
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


def parse_vcf(vcf, args, name=None):
    """
    Parse a VCF file and return a list of sequences
    :param vcf: path to the VCF file
    :param name: name of the sample
    """
    calls = []
    if name is None:
        name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        coords = f"{v.CHROM}:{v.POS}"
        # args.repeat, if set, specifies coordinates of a repeat to plot
        if args.repeat and coords != args.repeat:
            continue
        sequences = parse_alts(v.ALT, v.genotypes[0])
        if args.somatic:
            if v.INFO.get("SEQS") is None:
                sys.stderr.write("WARNING: No SEQS field found in VCF, skipping\n")
                return []
            else:
                somatic_sequences = v.INFO.get("SEQS").split(",")
        if sequences[0] and len(sequences[0]) > args.minlen:
            if args.somatic:
                if len(somatic_sequences) > 0:
                    for i, s in enumerate(somatic_sequences[0].split(":")):
                        if len(s) > args.minlen:
                            calls.append((coords, name, f"Allele1_{i}", s))
                # if there are outliers, add them as well, which is done only once as those are not phased
                if v.INFO.get("OUTLIERS") is not None:
                    for i, s in enumerate(v.INFO.get("OUTLIERS").split(",")):
                        if len(s) > args.minlen:
                            calls.append((coords, name, f"Outlier_{i}", s))
            else:
                calls.append((coords, name, "Allele1", sequences[0]))
        if sequences[1] and len(sequences[1]) > args.minlen:
            if args.somatic:
                if len(somatic_sequences) > 1:
                    for i, s in enumerate(somatic_sequences[1].split(":")):
                        if len(s) > args.minlen:
                            calls.append((coords, name, f"Allele2_{i}", s))
            else:
                calls.append((coords, name, "Allele2", sequences[1]))
    if args.longest_only:
        if args.somatic:
            sys.exit("ERROR: --longest_only is not supported with --somatic")
        # if two alleles pass the minimal length cutoff, only keep the longest one
        if calls:
            calls = [max(calls, key=lambda x: len(x[3]))]
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
    parser = ArgumentParser(
        description="Create a repeat sequence plot similar to the pathSTR <sequence> composition vizualization, but stand-alone"
    )
    parser.add_argument(
        "--names", help="Sample names to use, comma-separated", default=None
    )
    parser.add_argument(
        "-k",
        "--kmer",
        help="Kmer length to use for plot, ignored with --motifs",
        default=3,
        type=int,
    )
    parser.add_argument(
        "-n",
        "--number",
        help="Number of kmers to plot, ignored with --motifs",
        default=10,
        type=int,
    )
    parser.add_argument(
        "--motifs",
        help="Manually specify the motifs to plot, comma separated",
        default=None,
    )
    parser.add_argument(
        "--repeat",
        help="Chromosome and POS of repeat to plot from VCF e.g. chr1:12345, default: all repeats",
        default=None,
    )
    parser.add_argument(
        "-o", "--out", help="Output file name (html)", default="astronaut.html"
    )
    parser.add_argument(
        "-m", "--minlen", help="Minimal allele length to plot", default=20, type=int
    )
    parser.add_argument("--hide-labels", help="Hide sample labels", action="store_true")
    parser.add_argument(
        "--label_size", help="Size of sample labels", default=8, type=int
    )
    parser.add_argument(
        "--alphabetic", help="Sort samples alphabetically", action="store_true"
    )
    parser.add_argument(
        "--hide_allele_label", help="Hide 'Allele' from labels", action="store_true"
    )
    parser.add_argument(
        "--publication",
        help="Create a plot suitable for publication",
        action="store_true",
    )
    parser.add_argument(
        "--title", help="Title of the plot", default="Repeat composition"
    )
    parser.add_argument(
        "--sampleinfo", help="TSV file with sample information", default=None
    )
    parser.add_argument(
        "--somatic",
        help="Parse somatic VCFs, expects a SEQS format field",
        action="store_true",
    )
    parser.add_argument("--size", help="Size of markers to plot", default=3, type=int)
    parser.add_argument(
        "--longest_only",
        help="Only plot the longest allele per individual",
        action="store_true",
    )
    parser.add_argument(
        "--legend_corner",
        help="Corner of the legend",
        default="bottomright",
        choices=["topright", "bottomright"],
    )
    parser.add_argument("--height", help="Height of the plot", default=800, type=int)
    parser.add_argument(
        "-t",
        "--table",
        help="Use a tsv (with name, sequence and group columns) as input",
    )
    parser.add_argument("vcf", help="VCF files to analyze", nargs="*")
    args = parser.parse_args()
    if args.names:
        if len(args.names.split(",")) != len(args.vcf):
            sys.exit(
                "ERROR: Number of names does not match number of VCFs\nNames: {}\nVCFs: {}".format(
                    args.names.split(","), args.vcf
                )
            )
    if args.table and args.vcf:
        sys.exit("ERROR: Please provide either VCFs or a table, not both")
    if not args.table and not args.vcf:
        sys.exit("ERROR: Please provide either VCFs or a table")
    return args


if __name__ == "__main__":
    main()
