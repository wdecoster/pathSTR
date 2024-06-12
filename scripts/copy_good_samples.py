from argparse import ArgumentParser
import pandas as pd
import shutil
import os


def main():
    args = get_args()
    cramino = pd.read_csv(args.cramino, sep="\t")
    good_samples = (
        cramino.loc[cramino["Yield [Gb]"] > 32, "identifier"]
        .str.replace(".hg38", "", regex=False)
        .str.split("_")
        .str[0]
        .str.split("-")
        .str[0]
        .str.replace("GM", "NA", regex=False)
        .tolist()
    )
    # remove the suspected Klinefelter case
    if "HG02372" in good_samples:
        good_samples.remove("HG02372")
    for vcf in args.variants:
        if vcf.split("/")[-1].replace(".vcf.gz", "") in good_samples:
            base = os.path.basename(vcf)
            shutil.copy(vcf, os.path.join(args.outdir, base))
    with open(os.path.join(args.outdir, "good_samples.txt"), "w") as f:
        f.write("\n".join(good_samples))


def get_args():
    parser = ArgumentParser("")
    parser.add_argument("-c", "--cramino", help="cramino overview table")
    parser.add_argument("-v", "--variants", help="vcf files from strdust", nargs="+")
    parser.add_argument(
        "-o", "--outdir", help="output directory", default="good_samples"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
