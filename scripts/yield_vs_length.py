from argparse import ArgumentParser
import plotly.express as px
import pandas as pd
import plotly


def main():
    args = get_args()
    plot_yield_vs_length(args.input, args.output)


def plot_yield_vs_length(cramino_file, output_file):
    df = pd.read_csv(cramino_file, sep="\t")
    fig = px.scatter(
        df,
        x="N50",
        y="Yield [Gb]",
        hover_data=["identifier"],
        marginal_y="violin",
    )
    fig.add_hline(y=32, line_dash="dot", line_color="red")
    fig.update_xaxes(ticks="inside", dtick=10000)
    fig.update_traces(marker=dict(size=3))

    fig.update_layout(
        xaxis_title="N50 [bp]",
        yaxis_title="Yield [Gb]",
        plot_bgcolor="white",
        width=500,
        height=800,
    )
    fig = fix_marginal_violin(fig)
    with open(output_file, "w") as f:
        f.write(fig.to_html())


def fix_marginal_violin(fig):
    for trace in fig.data:
        if type(trace) == plotly.graph_objs._violin.Violin:
            trace.points = False
            trace.spanmode = "hard"
    return fig


def get_args():
    parser = ArgumentParser("")
    parser.add_argument("-i", "--input", help="Merged cramino file")
    parser.add_argument(
        "-o", "--output", help="Output plot file", default="yield_vs_length.html"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
