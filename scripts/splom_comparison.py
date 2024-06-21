import gzip
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
from os import path
import numpy as np
import plotly.graph_objects as go


def main():
    # Get arguments
    args = get_args()
    # Make dictionary to store repeat values
    res = []
    # Get sample IDs from STRdust VCF file names
    names = [path.basename(v).replace(".vcf.gz", "") for v in args.strdust]
    # Open STRchive bed file
    columns = ["chrom", "begin"]
    df = pd.read_csv(args.bed, sep="\t", usecols=[0, 1], names=columns)
    # Open VCF files from STRdust and LongTR
    for caller, tag in zip(
        [args.strdust, args.longtr, args.benchmark], ["RB", "GB", "SVLEN"]
    ):
        genotypes = defaultdict(list)
        if caller != None:
            for vcf_file, name in zip(caller, names):
                f = (
                    gzip.open(vcf_file, "rt")
                    if vcf_file.endswith(".gz")
                    else open(vcf_file, "r")
                )
                for line in f:
                    if line.startswith("#"):
                        continue  # Skip header lines
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    pos = fields[1]
                    # Extract allele length values, either from INFO or FORMAT field
                    if tag == "SVLEN":
                        info = fields[7]
                        repeat_values = extract_svlen(info, tag)
                    elif tag == "RB" or tag == "GB":
                        format = fields[8]
                        sample_str = fields[9]
                        ref = fields[3]
                        repeat_values = extract_rb_gb(format, sample_str, tag, ref)
                    # Make sure . is replaced by None
                    repeat_values = [None if i == "." else i for i in repeat_values]
                    # Convert repeat values from str to float and sort them
                    if len(repeat_values) == 2:
                        repeat_values = [
                            float(i) if i != None else None for i in repeat_values
                        ]
                        if repeat_values[0] != None and repeat_values[1] != None:
                            repeat_values.sort()
                    elif len(repeat_values) == 1:
                        if repeat_values[0] != None:
                            repeat_values = [float(repeat_values[0]), None]
                        else:
                            repeat_values = [None, None]
                    else:
                        repeat_values = [None, None]
                    # Find closest position in STRchive bed file
                    pos_bed = find_closest_position(df, chrom, int(pos))
                    # Store repeat values in dictionary
                    # Done this way because benchmark VCF file sometimes has multiple lines for the same locus
                    locus = f"{chrom}:{pos_bed}_{name}"
                    # Check if key_1 already contains a value
                    if len(genotypes[f"{locus}_1"]) == 0:
                        genotypes[f"{locus}_1"].append(repeat_values[0])
                        # If there is a second value, add it to the second allele
                        # Benchmark VCF file only contains 1 value per line
                        if len(repeat_values) > 1 and tag != "SVLEN":
                            genotypes[f"{locus}_2"].append(repeat_values[1])
                    # If key_1 already contains a value, store the repeat value in key_2
                    else:
                        genotypes[f"{locus}_2"].append(repeat_values[0])
                    # Make sure that the repeat values are stored in the right order (smallest first, None last)
                    if (
                        len(genotypes[f"{locus}_1"]) > 0
                        and len(genotypes[f"{locus}_2"]) > 0
                    ):
                        if (
                            genotypes[f"{locus}_1"][0] == None
                            or genotypes[f"{locus}_2"][0] == None
                        ):
                            if (
                                genotypes[f"{locus}_1"][0] == None
                                and genotypes[f"{locus}_2"][0] != None
                            ):
                                genotypes[f"{locus}_1"][0] = genotypes[f"{locus}_2"][0]
                                genotypes[f"{locus}_2"][0] = None
                        elif genotypes[f"{locus}_1"][0] > genotypes[f"{locus}_2"][0]:
                            genotypes[f"{locus}_1"][0], genotypes[f"{locus}_2"][0] = (
                                genotypes[f"{locus}_2"][0],
                                genotypes[f"{locus}_1"][0],
                            )
                f.close()
                # Add missing values (None) to dictionary, so that all arrays have the same length
                for key in genotypes:
                    if len(genotypes[key]) == 0:
                        genotypes[key].append(None)
            # Make dataframe
            df_output = pd.DataFrame.from_dict(genotypes)
            # Transpose the dataframe, rename columns, and set the chrom:pos_sample column as the index
            df_pivoted = df_output.transpose().reset_index()
            if tag == "RB":
                df_pivoted.columns = ["chrom:pos_sample", "STRdust"]
            elif tag == "GB":
                df_pivoted.columns = ["chrom:pos_sample", "LongTR"]
            elif tag == "SVLEN":
                df_pivoted.columns = ["chrom:pos_sample", "Benchmark"]
            df_pivoted.set_index("chrom:pos_sample", inplace=True)
            res.append(df_pivoted)

    # Merge dataframes
    df_merged = res[0]
    for i in range(1, len(res)):
        df_merged = pd.merge(
            df_merged, res[i], how="outer", left_index=True, right_index=True
        )

    # Make sure that the benchmark values are stored at the right allele
    if "Benchmark" in df_merged.columns:
        for index, row in df_merged.iterrows():
            if index.endswith("_1"):
                locus_2 = index.replace("_1", "_2")
                row_2 = df_merged.loc[locus_2]
                # If both alleles have a benchmark value, they are already sorted from small to large
                if row["Benchmark"] != None and row_2["Benchmark"] != None:
                    continue
                # If only one allele has a benchmark value, store it in the right allele
                elif row["Benchmark"] != None:
                    diff_1 = []
                    diff_2 = []
                    # Gather absolute differences to allele 1 and 2
                    for caller, tool in zip(
                        [args.strdust, args.longtr], ["STRdust", "LongTR"]
                    ):
                        if caller != None:
                            if row[tool] != None:
                                diff_1.append(abs(row["Benchmark"] - row[tool]))
                            if row_2[tool] != None:
                                diff_2.append(abs(row["Benchmark"] - row_2[tool]))
                    # Calculate mean difference to allele 1 and 2
                    mean_1 = np.mean(diff_1)
                    mean_2 = np.mean(diff_2)
                    # If the mean difference to allele 2 is smaller, switch the values
                    if mean_2 < mean_1:
                        df_merged.at[locus_2, "Benchmark"] = row["Benchmark"]
                        df_merged.at[index, "Benchmark"] = None

    # Calculate and save R² values
    r_squared_results = calculate_r_squared(df_merged)

    # Make scatterplot matrix
    splom = make_splom(df_merged, r_squared_results, args)

    # Save scatterplot matrix
    splom.write_html(args.output)


def get_args():
    parser = ArgumentParser(description="scatterplot matrix for tandem repeats")
    parser.add_argument("--strdust", help="VCF(s) from STRdust", nargs="+")
    parser.add_argument("--longtr", help="VCF(s) from LongTR", nargs="*")
    parser.add_argument("--benchmark", help="HG002 benchmark VCF", nargs="*")
    parser.add_argument("--bed", help="STRchive bed file", required=True)
    parser.add_argument("--output", help="Output html file name", required=True)
    args = parser.parse_args()
    return args


# Function to parse INFO field and extract SVLEN values
def extract_svlen(info, tag):
    # Find the tag in the INFO field
    info_fields = info.split(";")
    tag_field = [field for field in info_fields if field.startswith(f"{tag}=")]
    # If the tag is found, extract the repeat values
    if len(tag_field) == 1:
        repeat_values = tag_field[0].split("=")[1]
        # If the tag is SVLEN, check if the SVTYPE is DEL (-> multiply the repeat values by -1)
        if tag == "SVLEN":
            svtype_field = [
                field for field in info_fields if field.startswith("SVTYPE=")
            ]
            if len(svtype_field) == 1:
                svtype = svtype_field[0].split("=")[1]
                if svtype == "DEL":
                    repeat_values = float(repeat_values) * -1
            repeat_values = [repeat_values]
        return repeat_values
    else:
        return [None, None]


# Function to extract RB or GB values
def extract_rb_gb(format, sample_str, tag, ref):
    # Find at which location the tag is in the format field
    format_fields = format.split(":")
    i = format_fields.index(tag)
    # Get the corresponding values from the sample string
    sample_fields = sample_str.split(":")
    al_gb_field = sample_fields[i]
    # RB values are formatted as x,y and GB values are formatted as x|y
    if tag == "RB":
        repeat_values = al_gb_field.split(",")
    elif tag == "GB":
        repeat_values = al_gb_field.split("|")
    return repeat_values


# Function to find the closest position in the original STRchive bed file
def find_closest_position(df, input_chrom, input_pos):
    # Filter the DataFrame for the input chromosome
    df_chrom = df[df["chrom"] == input_chrom].reset_index(drop=True)
    # If there is only one STR position on the chromosome, return that position
    if len(df_chrom) == 1:
        return df_chrom["begin"].iloc[0]
    # Calculate the absolute difference between input position and start positions
    df_chrom["difference"] = abs(df_chrom["begin"] - input_pos)
    # Find the row with the smallest difference
    closest_row = df_chrom.loc[df_chrom["difference"] == df_chrom["difference"].min()]
    # Warning if the closest position is more than 1000 bp away
    if closest_row["difference"].iloc[0] > 1000:
        print(
            f"Warning: {input_chrom}:{input_pos} is more than 1000 bp away from the closest STR position"
        )
    # Return the closest position of the bed file
    return closest_row["begin"].iloc[0]


# Function to calculate R² values between all pairs of software tools
def calculate_r_squared(df):
    columns = df.columns
    r_squared_values = {}
    for i in range(len(columns)):
        for j in range(i + 1, len(columns)):
            tool1 = columns[i]
            tool2 = columns[j]
            # Filter out rows with missing values in either tool1 or tool2
            valid_rows = (
                df[[tool1, tool2]].apply(pd.to_numeric, errors="coerce").dropna()
            )
            # Calculate the correlation matrix and extract the R² value
            correlation_matrix = np.corrcoef(valid_rows[tool1], valid_rows[tool2])
            r_squared = correlation_matrix[0, 1] ** 2
            # Store the R² value in a dictionary
            r_squared_values[f"{tool1} vs {tool2}"] = r_squared
    # Create a dataframe from the dictionary
    r_squared_df = pd.DataFrame(
        list(r_squared_values.items()), columns=["Software tool comparison", "R² Value"]
    )
    return r_squared_df


def make_splom(df_merged, r_squared_results, args):
    # Make dimensions for scatter plot matrix: only add tools that had input files
    dimensions = []
    for caller, tool in zip(
        [args.strdust, args.longtr, args.benchmark], ["STRdust", "LongTR", "Benchmark"]
    ):
        if caller != None:
            dimensions.append(dict(label=tool, values=df_merged[tool]))
    # Make scatter plot matrix
    fig = go.Figure(
        go.Splom(
            dimensions=dimensions,
            diagonal_visible=False,
            showupperhalf=False,
            marker=dict(color="black", size=5),
            hoverinfo="text",
            text=df_merged.index,
        )
    )
    fig.update_layout(
        plot_bgcolor="#eeeeee",
        yaxis=dict(showgrid=False),
        yaxis2=dict(showgrid=False),
        xaxis=dict(showgrid=False),
        xaxis2=dict(showgrid=False),
    )
    # Find highest value of repeat length per tool to be able to locate the R² annotations
    max_values = {}
    for dim in dimensions:
        max_value = df_merged[dim["label"]].max()
        max_values[dim["label"]] = max_value
    # Add annotations for R² values
    for i, dim1 in enumerate(fig.data[0]["dimensions"]):
        for j, dim2 in enumerate(fig.data[0]["dimensions"]):
            if i != j:
                col1 = dim1["label"]
                col2 = dim2["label"]
                name = f"{col1} vs {col2}"
                r_squared = r_squared_results.loc[
                    r_squared_results["Software tool comparison"] == name, "R² Value"
                ].values
                if r_squared.size != 0:
                    fig.add_annotation(
                        # Locate R² annotations relatively to the lengths of the axes
                        x=float(max_values[col1]) * 0.6,
                        y=float(max_values[col2]) * 0.95,
                        text=f"R²={r_squared[0]:.3f}",
                        showarrow=False,
                        font=dict(size=12),
                        xref="x{}".format(i + 1),
                        yref="y{}".format(j + 1),
                    )
    # Add a background box to the R² annotations
    for i in range(len(fig.layout.annotations)):
        fig.layout.annotations[i].update(
            bgcolor="rgba(255,255,255,0.7)", font=dict(color="black")
        )
    # Return figure
    return fig


if __name__ == "__main__":
    main()
