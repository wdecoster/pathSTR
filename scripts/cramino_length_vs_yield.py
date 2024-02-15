import pandas as pd
import plotly.express as px
from argparse import ArgumentParser


def main():
    args = get_args()
    res = []
    for cramino_file in args.input:
        for line in open(cramino_file):
            name = cramino_file.split("/")[-1].split(".")[0]
            if line.startswith("Yield"):
                yield_gb = float(line.rstrip().split()[-1])
            if line.startswith("N50"):
                n50 = int(line.rstrip().split()[-1])
            res.append((name, yield_gb, n50))
    df = pd.DataFrame(res, columns=["sample", "yield", "n50"])
    fig = px.scatter(df, x="n50", y="yield", title="cramino N50 vs yield [Gb]")
    fig.update_layout(
        xaxis_title="N50",
        yaxis_title="Yield [Gb]",
    )
    print(fig.to_html())


def get_args():
    parser = ArgumentParser(
        description="Plot the length of the repeats against the yield of the reads"
    )
    parser.add_argument("-i", "--input", help="Cramino files", nargs="+", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    main()
