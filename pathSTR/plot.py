import sys
import plotly.express as px
from pathSTR.count_kmers import parse_kmers
from plotly.subplots import make_subplots
from math import ceil, log10
import plotly.graph_objects as go
import plotly


def violin_plot(
    filtered_df, repeats, selected_gene, violin_options=None, publication_ready=False
):
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
        title_text=f"Length distribution of {selected_gene} repeat",
    )
    if filtered_df["Group"].nunique() > 1 and "sex" not in violin_options:
        fig.update_layout(legend_title_text="Group")
    elif "sex" in violin_options:
        fig.update_layout(legend_title_text="Sex")
    else:
        fig.update_layout(showlegend=False)
    if "pathlen" in violin_options:
        # path length has to be corrected for the reference length
        # pathogenic length is in motif units, so the reflen has to be converted to motif units too
        path_length = (
            repeats.pathogenic_min_length(selected_gene)
            - repeats.reflen(selected_gene, unit="units")
            if "ref_diff" in violin_options
            else repeats.pathogenic_min_length(selected_gene)
        )
        # if path_length is larger than the current y-axis, extend the y-axis
        if path_length > filtered_df["length"].max():
            if "log" in violin_options:
                fig.update_layout(yaxis_range=[0, log10(ceil(path_length * 1.1))])
            else:
                fig.update_layout(yaxis_range=[1, ceil(path_length * 1.1)])
        fig.add_hline(y=path_length, line_dash="dot", line_color="red")
    if publication_ready:
        fig.update_layout(
            font=dict(size=20),
            legend=dict(
                title_font=dict(size=16),
                font=dict(size=16),
            ),
            plot_bgcolor="white",
            width=800,
            height=800,
        )
        fig.update_traces(marker=dict(size=8, opacity=0.8))
    return fig


def length_scatter(
    filtered_df,
    selected_gene,
    path_length=None,
    violin_options=None,
    publication_ready=False,
):
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
        .astype(  # cast object to float to avoid futurewarning with fillna causing a downcast
            {
                "length_Allele1": "float64",
                "length_Allele2": "float64",
                "ref_diff_Allele1": "float64",
                "ref_diff_Allele2": "float64",
            }
        )
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
        title=f"Lengths per allele of {selected_gene} repeat",
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

    if publication_ready:
        fig.update_layout(
            font=dict(size=20),
            legend=dict(
                title_font=dict(size=16),
                font=dict(size=16),
            ),
            plot_bgcolor="white",
            width=800,
            height=800,
        )
        for trace in fig.data:
            if type(trace) == plotly.graph_objs.scatter:
                trace.marker.update(dict(size=8, opacity=0.8))
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


def kmer_plot_raw(
    kmer_df,
    selected_gene,
    length_range=None,
    sort=None,
    num_columns=0,
    publication_ready=False,
):
    """
    Create a heatmap of the kmer counts per allele and per sample

    """
    if length_range:
        min_length = length_range[0]
        max_length = length_range[1]
        kmer_df = kmer_df[
            kmer_df["length"].between(min_length, max_length, inclusive="both")
        ].drop(columns=["length"])
    else:
        kmer_df = kmer_df.drop(columns=["length"])
    # sorting the columns based on their total frequency
    kmer_frequency_sorted = (
        kmer_df.sum(axis="index").sort_values(ascending=False).index.to_list()
    )
    kmer_df = kmer_df[kmer_frequency_sorted]

    if sort:
        # in the raw mode, the kmer_options indicates the kmer(s) to sort the rows on
        kmer_df = kmer_df.sort_values(by=sort, ascending=False)
    min_seen_kmer = 0.98 * len(kmer_df.index)
    # only keep columns that are not < 0.01 for too many samples
    mask1 = (kmer_df < 0.01).sum(axis=0) < (min_seen_kmer)
    # but keep also columns that are above 0.1 for at least one sample
    mask2 = (kmer_df > 0.1).sum(axis=0) > 0

    kmer_df = kmer_df.loc[:, mask1 | mask2]

    # the plot takes up a terrible lot of vertical space, so try splitting it up in <columns> columns
    # this however could be a problem on a small screen
    # maybe add another slider to select the number of columns
    if num_columns:
        columns = num_columns
    else:
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
                color_continuous_scale=[(0, "white"), (1, "darkblue")],
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
            title_text=f"Kmer frequency for {selected_gene} repeat",
            plot_bgcolor="rgba(0, 0, 0, 0)",
            paper_bgcolor="rgba(0, 0, 0, 0)",
            height=height,
            width=height / 4,
        )
    )
    fig.update_coloraxes(colorscale="Blues")
    if publication_ready:
        fig.update_layout(
            font=dict(size=16),
            legend=dict(
                title_font=dict(size=16),
                font=dict(size=16),
            ),
            plot_bgcolor="white",
            width=800,
            height=800,
        )
    return fig


def kmer_plot_collapsed(
    kmer_df,
    selected_gene,
    length_range=None,
    min_group_size=0,
    publication_ready=False,
):
    if length_range:
        min_length = length_range[0]
        max_length = length_range[1]
        kmer_df = kmer_df[
            kmer_df["length"].between(min_length, max_length, inclusive="both")
        ].drop(columns=["length"])
    else:
        kmer_df = kmer_df.drop(columns=["length"])
    # sorting the columns based on their total frequency
    kmer_frequency_sorted = (
        kmer_df.sum(axis="index").sort_values(ascending=False).index.to_list()
    )
    kmer_df = kmer_df[kmer_frequency_sorted]
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
    collapsed = collapsed[collapsed["count"] > min_group_size]

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
            "title_text": f"Kmer frequency for {selected_gene} repeat",
            "plot_bgcolor": "white",
            "paper_bgcolor": "white",
            "height": 1200,
            "width": 1200,
        }
    )

    fig["layout"]["xaxis2"].update(title="Number of carriers")
    fig.update_coloraxes(colorscale=[(0, "white"), (1, "darkblue")])
    if publication_ready:
        fig.update_layout(
            font=dict(size=16),
            legend=dict(
                title_font=dict(size=16),
                font=dict(size=16),
            ),
            plot_bgcolor="white",
            width=800,
            height=800,
        )
        # increase font size of axis labels
        fig.update_xaxes(tickfont_size=12)
    return fig


def kmer_plot_sequence(
    kmer_df,
    repeat_df,
    selected_gene,
    pathogenic_length=None,
    length_range=None,
    direction="left-to-right",
    publication_ready=False,
):
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
    if length_range:
        min_length = length_range[0]
        max_length = length_range[1]
        kmer_df = kmer_df[
            kmer_df["length"].between(min_length, max_length, inclusive="both")
        ].drop(columns=["length"])
    else:
        kmer_df = kmer_df.drop(columns=["length"])
    # sorting the columns based on their total frequency
    kmer_frequency_sorted = (
        kmer_df.sum(axis="index").sort_values(ascending=False).index.to_list()
    )
    kmer_df = kmer_df[kmer_frequency_sorted]

    # using only the 10 most frequent kmers by limiting kmer_df to the first 10 columns
    # assign a color to each kmer
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

    # for the alt sequence of every individual,
    # plot the order of the 10 most frequent kmers, with others in grey
    if length_range:
        min_length, max_length = length_range
        repeat_df = (
            repeat_df[
                repeat_df["length"].between(min_length, max_length, inclusive="both")
            ]
            .sort_values(by="length", ascending=False)
            .dropna(subset=["sequence"])
        )
    else:
        repeat_df = repeat_df.sort_values(by="length", ascending=False).dropna(
            subset=["sequence"]
        )
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
    if direction == "right-to-left":
        # in the sequence mode, the kmer_options indicates the alignment of the plot relative to the x-axis
        # split the seq_colored into a list and reverse the list of colors
        repeat_df["seq_colored"] = (
            repeat_df["seq_colored"].str.split(";").apply(lambda x: x[::-1])
        )
    else:
        # split the seq_colored into a list
        repeat_df["seq_colored"] = repeat_df["seq_colored"].str.split(";")
    # add a range column to the dataframe, to enumerate the nucleotides
    repeat_df["range"] = repeat_df["seq_colored"].apply(lambda x: list(range(len(x))))
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
        category_orders={"identifier": repeat_df["identifier"][::-1]},
        title=f"Sequence of {selected_gene} repeat",
    )
    fig.update_traces(marker=dict(size=3))
    # make the y axis labels a bit smaller
    fig.update_yaxes(tickfont_size=8)
    if pathogenic_length:
        fig.add_vline(x=pathogenic_length, line_dash="dot", line_color="red")
    if direction == "right-to-left":
        fig.update_layout(xaxis_autorange="reversed")
    # set the height of the plot depending on the number of samples
    fig.update_layout(
        height=get_height(len(repeat_df)),
        width=get_width(repeat_df["sequence"].apply(lambda x: len(x)).max()),
    )
    if publication_ready:
        fig.update_layout(
            font=dict(size=16),
            legend=dict(
                title_font=dict(size=16),
                font=dict(size=16),
            ),
        )
    # change the size of the dots in the legend
    fig.update_layout(
        legend={"itemsizing": "constant"},
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )
    return fig


def get_width(longest_allele):
    """
    Return the width of the plot based on the length of the longest allele
    The width is constrained between 800 and 1600, and dynamic between those values
    """
    max_length = longest_allele * 10
    if max_length > 1600:
        return 1600
    elif max_length < 800:
        return 800
    else:
        return max_length


def get_height(num_samples):
    """
    Return the height of the plot based on the number of samples
    """
    return max(50, num_samples)
