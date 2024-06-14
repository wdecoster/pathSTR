# this script will take a directory of VCF files, as well as a genome reference build and genotyper, to output a file of filenames to stdout

from argparse import ArgumentParser
import glob


def main():
    args = get_args()
    for vcf in glob.glob(f"{args.vcf_dir}/*.vcf.gz"):
        print(f"{vcf}\t{args.build}\t{args.caller}")


def get_args():
    parser = ArgumentParser(description="Prepare a file of filenames for pathSTR")
    parser.add_argument(
        "--vcf_dir",
        help="Directory of VCF files",
    )
    parser.add_argument(
        "--build",
        help="Reference genome build",
    )
    parser.add_argument(
        "--caller",
        help="Genotyper",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
