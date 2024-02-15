import pandas as pd
from argparse import ArgumentParser


def main():
    args = get_args()
    cramino = pd.concat(
        [pd.read_csv(f, sep="\t").set_index("File name") for f in args.input],
        axis="columns",
    )
    cramino.columns = [
        c.replace(".cram", "").replace("_Prom", "") for c in cramino.columns
    ]
    cramino = cramino.transpose().reset_index(names="identifier")
    cramino.to_csv(args.output, sep="\t", index=False)


def get_args():
    parser = ArgumentParser("")
    parser.add_argument("-i", "--input", nargs="+", help="input files")
    parser.add_argument("-o", "--output", help="output files")
    return parser.parse_args()


if __name__ == "__main__":
    main()
