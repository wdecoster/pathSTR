import plotly.express as px


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
