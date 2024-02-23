import sys
import plotly.express as px
from pathSTR.count_kmers import parse_kmers
from plotly.subplots import make_subplots
from math import ceil, log10
import plotly.graph_objects as go


def violin_plot(filtered_df, path_length=None, violin_options=None):
    if "density" in violin_options:
        fig = px.violin(
            filtered_df,
            x="Superpopulation" if "population" in violin_options else "gene",
            y="ref_diff" if "ref_diff" in violin_options else "length",
            color="Sex" if "sex" in violin_options else "Group",
            log_y="log" in violin_options,
            points="all",
            hover_data=["sample"],
        )
        fig.update_traces(spanmode="hard", marker=dict(size=3))
    else:
        fig = px.strip(
            filtered_df,
            x="Superpopulation" if "population" in violin_options else "gene",
            y="ref_diff" if "ref_diff" in violin_options else "length",
            color="Sex" if "sex" in violin_options else "Group",
            log_y="log" in violin_options,
            hover_data=["sample"],
        )
        fig.update_traces(marker=dict(size=3))
    fig.update_layout(
        xaxis_title="",
        yaxis_title=(
            "Repeat length [log(units)]"
            if "log" in violin_options
            else "Repeat length [units]"
        ),
    )
    if filtered_df["Group"].nunique() > 1 and "sex" not in violin_options:
        fig.update_layout(legend_title_text="Group")
    elif "sex" in violin_options:
        fig.update_layout(legend_title_text="Sex")
    else:
        fig.update_layout(showlegend=False)
    if "pathlen" in violin_options:
        # if path_length is larger than the current y-axis, extend the y-axis
        if path_length > filtered_df["length"].max():
            if "log" in violin_options:
                fig.update_layout(yaxis_range=[0, log10(ceil(path_length * 1.1))])
            else:
                fig.update_layout(yaxis_range=[1, ceil(path_length * 1.1)])
        fig.add_hline(y=path_length, line_dash="dot", line_color="red")
    return fig


def length_scatter(filtered_df, path_length=None, violin_options=None):
    pivot_df = filtered_df.pivot(
        index="sample",
        columns="allele",
        values=["length", "ref_diff", "Sex", "Superpopulation", "Group"],
    )
    pivot_df.columns = [
        "_".join([str(v) for v in c]) for c in pivot_df.columns.to_flat_index()
    ]
    pivot_df = (
        pivot_df.drop(
            columns=["Sex_Allele2", "Superpopulation_Allele2", "Group_Allele2"]
        )
        .rename(
            columns={
                "Sex_Allele1": "Sex",
                "Superpopulation_Allele1": "Superpopulation",
                "Group_Allele1": "Group",
            }
        )
        .reset_index()
        .fillna(  # fill NaNs with 0 for the scatter plot of males on chrX
            {
                "ref_diff_Allele1": 0,
                "ref_diff_Allele2": 0,
                "length_Allele1": 0,
                "length_Allele2": 0,
            }
        )
    )
    pivot_df["longest_allele"] = pivot_df.apply(
        lambda x: max(x["length_Allele1"], x["length_Allele2"]), axis=1
    )
    pivot_df["shortest_allele"] = pivot_df.apply(
        lambda x: min(x["length_Allele1"], x["length_Allele2"]), axis=1
    )
    pivot_df["ref_diff_longest"] = pivot_df.apply(
        lambda x: max(x["ref_diff_Allele1"], x["ref_diff_Allele2"]), axis=1
    )
    pivot_df["ref_diff_shortest"] = pivot_df.apply(
        lambda x: min(x["ref_diff_Allele1"], x["ref_diff_Allele2"]), axis=1
    )

    fig = px.scatter(
        pivot_df,
        x="ref_diff_longest" if "ref_diff" in violin_options else "longest_allele",
        y="ref_diff_shortest" if "ref_diff" in violin_options else "shortest_allele",
        color="Sex" if "sex" in violin_options else "Group",
        symbol="Superpopulation" if "population" in violin_options else None,
        log_x="log" in violin_options,
        log_y="log" in violin_options,
        hover_data=[
            "sample",
            "longest_allele",
            "shortest_allele",
            "ref_diff_longest",
            "ref_diff_shortest",
        ],
    )
    fig.update_traces(marker=dict(size=3))
    if "density" in violin_options:
        fig.add_histogram2dcontour(
            x=(
                pivot_df["ref_diff_longest"]
                if "ref_diff" in violin_options
                else pivot_df["longest_allele"]
            ),
            y=(
                pivot_df["ref_diff_shortest"]
                if "ref_diff" in violin_options
                else pivot_df["shortest_allele"]
            ),
            xaxis="x",
            yaxis="y",
            colorscale="Blues",
            showscale=False,
        )
    if "pathlen" in violin_options:
        # if path_length is larger than the current y-axis, extend the y-axis
        if path_length > filtered_df["length"].max():
            if "log" in violin_options:
                fig.update_layout(yaxis_range=[0, log10(ceil(path_length * 1.1))])
                fig.update_layout(xaxis_range=[0, log10(ceil(path_length * 1.1))])
            else:
                fig.update_layout(yaxis_range=[1, ceil(path_length * 1.1)])
                fig.update_layout(xaxis_range=[1, ceil(path_length * 1.1)])
        fig.add_hline(y=path_length, line_dash="dot", line_color="red")
        fig.add_vline(x=path_length, line_dash="dot", line_color="red")
    else:
        if "log" in violin_options:
            fig.update_layout(
                xaxis_range=[0, log10(ceil(filtered_df["length"].max() * 1.1))]
            )
            fig.update_layout(
                yaxis_range=[0, log10(ceil(filtered_df["length"].max() * 1.1))]
            )
    if (
        pivot_df["Group"].nunique() > 1
        or "sex" in violin_options
        or "population" in violin_options
    ):
        fig.update_layout(legend_title_text="Group")
    else:
        fig.update_layout(showlegend=False)

    # make the y-range the same as the x-range

    fig.update_layout(
        height=1000,
        width=1000,
        xaxis_title="Repeat length longer allele [units]",
        yaxis_title="Repeat length shorter allele [units]",
        # use a white background
        plot_bgcolor="rgba(0, 0, 0, 0)",
    )
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


def kmer_plot(kmer_df, repeat_df, mode="collapsed", min_length=0, sort_by=None):
    """
    Create plots of kmers found in the repeat sequences.
    The mode can be "raw", "collapsed" or "sequence"
    raw: plot the heatmap of per repeat and per allele kmers - but that is a whole lot of data
    collapsed: group similar samples together and plot the heatmap with a marginal histogram
    sequences: plot the sequence of the 10 most frequent kmers in the order that they're found
    optionally the minimum expansion length can be set to filter out alleles of short repeats

    :param kmer_df: the dataframe with the kmer counts
    :param mode: the mode of the plot, either "raw", "collapsed" or "sequence"
    :param min_length: the minimum length of the repeat to be included in the plot
    :param sort_by: the kmer to sort the rows on (only for mode="raw")
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
        if sort_by:
            kmer_df = kmer_df.sort_values(by=sort_by, ascending=False)
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
                    kmer_df[i * batch_num : (i + 1) * batch_num][::-1],
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
        collapsed["count"] = collapsed["identifier"].apply(lambda x: len(x.split(",")))
        collapsed = collapsed.sort_values(by="count", ascending=False).reset_index(
            drop=True
        )
        # create a heatmap with a marginal bar chart to indicate the number of carriers
        fig = make_subplots(
            rows=1,
            cols=2,
            column_widths=[0.8, 0.2],
            shared_yaxes=True,
            subplot_titles=("Kmer frequency", "Number of carriers"),
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
                x=collapsed["count"][::-1],
                y=collapsed.index[::-1],
                orientation="h",
            ).data[0],
            row=1,
            col=2,
        )

        # make the axis labels a bit smaller
        fig.update_xaxes(tickfont_size=8, tickangle=45)
        fig.update_yaxes(visible=False, showticklabels=False)
        fig.update_layout(
            {
                "plot_bgcolor": "rgba(0, 0, 0, 0)",
                "paper_bgcolor": "rgba(0, 0, 0, 0)",
                "height": 1200,
                "width": 1200,
            }
        )
        fig["layout"]["xaxis2"].update(title="Number of carriers in group")
        fig.update_coloraxes(colorscale="Blues")
        return fig
    elif mode == "sequence":
        # using only the 10 most frequent kmers by limiting kmer_df to the first 10 columns
        # assign a color to each kmer
        colors = px.colors.qualitative.Plotly
        kmer_dict = {k: c for k, c in zip(kmer_df.columns[:10], colors)}
        # for the alt sequence of every individual,
        # plot the order of the 10 most frequent kmers, with others in grey
        repeat_df = repeat_df[repeat_df["length"] >= min_length]
        # draw a scatter plot showing for each sample and allele the order of the kmers in the sequence
    else:
        sys.stderr.write("Invalid mode for kmer plot\n")
        return None
