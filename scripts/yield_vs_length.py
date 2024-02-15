from argparse import ArgumentParser


def main():
    args = get_args()
    plot_yield_vs_length(args.input, args.output)


def plot_yield_vs_length(cramino_file, output_file):
    import plotly.express as px
    import pandas as pd

    df = pd.read_csv(cramino_file, sep="\t")
    fig = px.scatter(df, x="N50", y="Yield [Gb]", hover_data=["identifier"])
    fig.update_layout(
        title="Yield vs N50",
        xaxis_title="N50 [bp]",
        yaxis_title="Yield [Gb]",
    )
    with open(output_file, "w") as f:
        f.write(fig.to_html())


def get_args():
    parser = ArgumentParser("")
    parser.add_argument("-i", "--input", help="Merged cramino file")
    parser.add_argument(
        "-o", "--output", help="Output plot file", default="yield_vs_length.html"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
