import pandas as pd
from argparse import ArgumentParser
import plotly.express as px
import numpy as np


def main():
    args = get_args()
    df = pd.read_csv(
        args.cramino, sep="\t", usecols=["identifier", "chrX", "chrY", "Yield [Gb]"]
    ).assign(
        identifier=lambda d: d["identifier"]
        .str.replace(".hg38", "", regex=False)
        .str.split("_")
        .str[0]
        .str.split("-")
        .str[0]
        .str.replace("GM", "NA", regex=False)
    )
    df = df[df["Yield [Gb]"] > args.minyield]
    # add a little bit of random noise to avoid overlapping points
    # the noise is normally distributed with a mean of 0 and a standard deviation of 0.02
    df["chrX"] = df["chrX"] + np.random.normal(0.0, 0.02, len(df))
    df["chrY"] = df["chrY"] + np.random.normal(0.0, 0.02, len(df))

    sampleinfo = pd.read_table(
        args.sampleinfo, usecols=["sample", "Sex", "source"]
    ).rename(columns={"sample": "identifier"})
    df = df.merge(sampleinfo, on="identifier", how="left").fillna("unknown")
    fig = px.scatter(
        df,
        x="chrX",
        y="chrY",
        color="Sex",
        hover_data=["identifier", "Yield [Gb]", "source"],
    )
    # change opacity to 0.6
    fig.update_traces(opacity=0.8, marker_size=3)
    # make the background white
    fig.update_layout(
        plot_bgcolor="white",
        scattermode="group",
        scattergap=0.75,
        width=500,
        height=500,
    )
    # update the axis labels
    fig.update_xaxes(
        title_text="chrX normalized coverage",
    )
    fig.update_yaxes(
        title_text="chrY normalized coverage",
    )
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=1,
            xanchor="right",
            x=0.97,
            orientation="h",
            bordercolor="lightgray",
            borderwidth=1,
        ),
        legend_title_text="",
    )

    # Add circle to highlight the weird sample
    # first find its coordinates
    x, y = (
        df.loc[(df["chrX"] > 0.8) & (df["chrY"] > 0.3), ["chrX", "chrY"]]
        .values[0]
        .tolist()
    )

    fig.add_shape(
        type="circle",
        xref="x",
        yref="y",
        x0=x - 0.02,
        y0=y - 0.02,
        x1=x + 0.02,
        y1=y + 0.02,
        line_color="lightgray",
    )
    print(fig.to_html())


def get_args():
    parser = ArgumentParser(
        "Plot copy numbers from multiple samples for a single region"
    )
    parser.add_argument("-c", "--cramino", help="summary file of cramino")
    parser.add_argument("--sampleinfo", help="excel file with sample information")
    parser.add_argument(
        "--minyield",
        help="minimum yield in gigabases to filter the results on first",
        default=32,
        type=int,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
