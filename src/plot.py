import plotly.express as px
from count_kmers import parse_kmers
from plotly.subplots import make_subplots
from math import ceil


def violin_plot(filtered_df, log=False, violin_options=None):
    fig = px.violin(
        filtered_df,
        x="Superpopulation" if "population" in violin_options else "gene",
        y="ref_diff" if "ref_diff" in violin_options else "length",
        color="Sex" if "sex" in violin_options else "Group",
        points="all",
        hover_data=["sample"],
    )
    if log:
        fig.update_layout(
            yaxis_type="log",
            xaxis_title="",
            yaxis_title="Repeat length [log(units)]",
        )
        # fig.update_layout(xaxis_range=[1, filtered_df["length"].max()])
    else:
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
    fig.update_traces(marker=dict(size=3))
    if filtered_df["Group"].nunique() > 1 and "sex" not in violin_options:
        fig.update_layout(legend_title_text="Group")
    elif "sex" in violin_options:
        fig.update_layout(legend_title_text="Sex")
    else:
        fig.update_layout(showlegend=False)
    return fig


def create_strip_plot(strip_df, log=False):
    fig = px.strip(
        strip_df,
        x="gene",
        y="length",
        color="Group",
        stripmode="overlay",
        hover_data=["sample"],
    )
    if log:
        fig.update_layout(yaxis_type="log")
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [log(units)]")
    else:
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
    fig.update_traces(marker=dict(size=2))
    if strip_df["Group"].nunique() > 1:
        fig.update_layout(legend_title_text="Group")
    else:
        fig.update_layout(showlegend=False)
    return fig


def kmer_plot(kmer_df, mode="collapsed", min_length=0):
    """
    Create plots of kmers found in the repeat sequences.
    The mode can be "raw", "collapsed" or "sequence"
    raw: plot the heatmap of per repeat and per allele kmers - but that is a whole lot of data
    collapsed: group similar samples together and plot the heatmap with a marginal histogram
    sequences: plot the sequence of the 10 most frequent kmers in the order that they're found
    optionally the minimum expansion length can be set to filter out alleles of short repeats
    """
    if min_length:
        kmer_df = kmer_df[kmer_df["length"] >= min_length].drop(columns=["length"])
    else:
        kmer_df = kmer_df.drop(columns=["length"])
    # sorting the columns based on their total frequency
    kmer_frequency_sorted = (
        kmer_df.sum(axis="index").sort_values(ascending=False).index.to_list()
    )
    kmer_df = kmer_df[kmer_frequency_sorted]
    if mode == "raw":
        min_seen_kmer = 0.98 * len(kmer_df.index)
        # only keep columns that are not < 0.01 for too many samples
        mask1 = (kmer_df < 0.01).sum(axis=0) < (min_seen_kmer)
        # but keep also columns that are above 0.1 for at least one sample
        mask2 = (kmer_df > 0.1).sum(axis=0) > 0

        kmer_df = kmer_df.loc[:, mask1 | mask2]

        # the plot takes up a terrible lot of vertical space, so try splitting it up in <columns> columns
        # this however could be a problem on a small screen
        # maybe add another slider to select the number of columns
        columns = ceil(len(kmer_df) / 500)
        fig = make_subplots(rows=1, cols=columns)
        batch_num = ceil(len(kmer_df) / columns)
        for i in range(0, columns):
            fig.add_trace(
                px.imshow(
                    kmer_df[i * batch_num : (i + 1) * batch_num],
                    labels=dict(
                        x="kmers",
                        y="identifier",
                        color="Frequency",
                    ),
                    color_continuous_scale="Blues",
                ).data[0],
                row=1,
                col=i + 1,
            )
        # make the axis labels a bit smaller
        fig.update_xaxes(tickfont_size=8, tickangle=45, side="top")
        fig.update_yaxes(tickfont_size=8)

        height = 6000 if columns > 1 else len(kmer_df) * 12

        fig.update_layout(
            dict(
                plot_bgcolor="rgba(0, 0, 0, 0)",
                paper_bgcolor="rgba(0, 0, 0, 0)",
                height=height,
                width=height / 4,
            )
        )
        fig.update_coloraxes(colorscale="Blues")
        return fig
    elif mode == "collapsed":
        # kmers are thrown out faster
        min_seen_kmer = 0.90 * len(kmer_df.index)
        # only keep columns that are not < 0.01 for too many samples
        mask1 = (kmer_df < 0.05).sum(axis=0) < (min_seen_kmer)
        # but keep also columns that are above 0.1 for at least 5 samples
        mask2 = (kmer_df > 0.1).sum(axis=0) > 5

        kmer_df = kmer_df.loc[:, mask1 | mask2]

        # round every kmer count to the nearest 0.1
        kmer_df = kmer_df.map(lambda x: round(x * 10) / 10)
        # group all samples together that have identical kmers numbers
        collapsed = (
            kmer_df.reset_index()
            .groupby(kmer_df.columns.to_list())["identifier"]
            .apply(",".join)
            .reset_index()
        )
        collapsed.to_csv("collapsed.csv", index=False, sep="\t")
        collapsed["count"] = collapsed["identifier"].apply(lambda x: len(x.split(",")))
        collapsed = collapsed.sort_values(by="count", ascending=False).reset_index(
            drop=True
        )
        # create a heatmap with a marginal bar chart to indicate the number of carriers
        fig = make_subplots(
            rows=1,
            cols=2,
            column_widths=[0.8, 0.2],
            subplot_titles=("Kmer frequency", "Frequency histogram"),
        )

        fig.add_trace(
            px.imshow(
                collapsed.drop(columns=["identifier", "count"]),
                labels=dict(x="kmers", y="group", color="Frequency"),
                color_continuous_scale="Blues",
            ).data[0],
            row=1,
            col=1,
        )
        fig.add_trace(
            px.bar(
                y=collapsed["count"],
                x=collapsed.index,
                orientation="h",
            ).data[0],
            row=1,
            col=2,
        )
        # make the axis labels a bit smaller
        fig.update_xaxes(tickfont_size=8, tickangle=45, side="top")
        fig.update_yaxes(tickfont_size=8)
        fig.update_layout(
            {
                "plot_bgcolor": "rgba(0, 0, 0, 0)",
                "paper_bgcolor": "rgba(0, 0, 0, 0)",
                "height": 800,
            }
        )
        fig.update_coloraxes(colorscale="Blues")
        return fig
