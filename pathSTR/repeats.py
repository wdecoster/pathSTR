import pandas as pd
import logging
import sys


class Repeats(object):
    def __init__(self, df=None, build=["hg38", "t2t"]):
        if df:
            self.df = df
        else:
            self.df = pd.concat(self.get_repeat_info(b) for b in build)

    def get_repeat_info(self, build):
        """
        This function parses a bed file as obtained from STRchive, which is intended for TRGT, but also works for STRdust.
        Also a CSV is downloaded, with extra information about the repeats.
        A dataframe is constructed, which is the data underlying the Repeats class instance

        The motif length that is used for e.g. kmer plots is the longest motif, which may have undesired consequences.
        But the shortest or random choice also was a bad thing, e.g. for FXN.
        """
        urls = {
            "hg38": "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/hg38.STRchive-disease-loci.TRGT.bed",
            "t2t": "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/T2T-chm13.STRchive-disease-loci.TRGT.bed",
        }
        bed = pd.read_csv(
            urls[build], sep="\t", header=None, names=["chrom", "start", "end", "info"]
        )
        bed["reflen"] = bed["end"] - bed["start"]
        # isolating the id to match up with the repeat info csv that is downloaded from STRchive
        bed["id"] = bed["info"].apply(
            lambda x: [i for i in x.split(";") if i.startswith("ID=")][0].replace(
                "ID=", ""
            )
        )
        # download the csv file from STRchive, and join it with the bed file to get the pathogenic_min length column
        bed = bed.join(
            pd.read_csv(
                "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/STRchive-database.csv",
                usecols=["id", "pathogenic_min"],
            ).set_index("id"),
            on="id",
        )

        # the TRGT bed file has as ID the <disease>_<gene> and we only care about the gene for pathSTR
        bed["name"] = bed["id"].apply(lambda x: x.split("_")[1])
        # in case there are duplicates in the name column, use the original ID from the info field
        dups = bed.duplicated(subset="name", keep=False)
        bed.loc[dups, "name"] = bed.loc[dups, "info"].apply(
            lambda x: [
                Repeats.fix_name(i) for i in x.split(";") if i.startswith("ID=")
            ][0].replace("ID=", "")
        )
        # extracting the motifs and their length based on the bed info field
        bed["motifs"] = bed["info"].apply(self.get_motifs)
        bed["motif_length"] = bed["motifs"].apply(self.get_longest_motif_length)
        # currently overwriting the id field of the table, but it would make more sense to have a more descriptive column name, e.g. "coords"
        bed["id"] = (
            bed["chrom"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str)
        )
        bed["build"] = build
        return bed.set_index("id")

    @staticmethod
    def fix_name(name):
        return f"{name.split('_')[1]}_{name.split('_')[0]}"

    @staticmethod
    def get_motifs(info):
        """
        For every info field from the bed file, get the value at the MOTIFS key
        """
        return [
            i.replace("MOTIFS=", "").split(",")
            for i in info.split(";")
            if i.startswith("MOTIFS=")
        ][0]

    @staticmethod
    def get_longest_motif_length(motifs):
        """
        Return the length of the longest motif
        """
        return sorted([len(m) for m in motifs], reverse=True)[0]

    def query(self, dataset, gene, column):
        """
        Generic method to query the dataframe based on a dataset, gene and column (or list of columns)
        """
        build = dataset.split("_")[1]
        try:
            return self.df.loc[
                (self.df["build"] == build) & (self.df["name"] == gene), column
            ].values[0]
        except IndexError:
            logging.error(f"Query in repeats failed for {build}, {gene}, {column}")
            logging.error(f"Dumping repeats object to repeats.tsv for debugging")
            self.df.to_csv("repeats.tsv", sep="\t")
            sys.exit(1)

    def motif_length(self, gene, dataset=None, build=None):
        if not build and not dataset:
            raise ValueError("Either build or dataset must be given")
        if dataset:
            return self.query(dataset, gene, "motif_length")
        else:
            # use dummy dataset name
            return self.query(f"dataset_{build}", gene, "motif_length")

    def motifs(self, gene, dataset):
        return self.query(dataset, gene, "motifs")

    def coords_to_gene(self, id, build):
        try:
            return self.df[self.df["build"] == build].loc[id, "name"]
        except KeyError:
            return None

    def gene_to_coords(self, gene, dataset):
        build = dataset.split("_")[1]
        return self.df.loc[
            (self.df["name"] == gene) & (self.df["build"] == build)
        ].index[0]

    def pathogenic_min_length(self, gene, dataset):
        return self.query(dataset, gene, "pathogenic_min")

    def reflen(self, gene, dataset, unit="units"):
        if unit == "units":
            return self.query(dataset, gene, "reflen")
        elif unit == "bp":
            return self.query(dataset, gene, "reflen") * self.motif_length(
                gene, dataset
            )
        else:
            raise ValueError("unit must be 'bp' or 'units'")

    def start(self, gene, dataset):
        return self.query(dataset, gene, "start")

    def end(self, gene, dataset):
        return self.query(dataset, gene, "end")

    def chrom(self, gene, dataset):
        return self.query(dataset, gene, "chrom")

    def coords(self, gene, dataset):
        return self.query(dataset, gene, ["chrom", "start", "end"])
