import pandas as pd


class Repeats(object):
    def __init__(self, bed=None, df=None):
        self.df = self.get_repeat_info()

    def get_repeat_info(self):
        """
        This function parses a bed file as obtained from STRchive, which is intended for TRGT, but also works for STRdust.
        Also a CSV is downloaded, with extra information about the repeats.
        A dataframe is constructed, which is the data underlying the Repeats class instance

        The motif length that is used for e.g. kmer plots is the longest motif, which may have undesired consequences.
        But the shortest or random choice also was a bad thing, e.g. for FXN.
        """
        url = "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/hg38.STRchive-disease-loci.TRGT.bed"
        bed = pd.read_csv(
            url, sep="\t", header=None, names=["chrom", "start", "end", "info"]
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

    def motif_length(self, gene):
        return self.df.loc[self.df["name"] == gene, "motif_length"].values[0]

    def motifs(self, gene):
        return self.df.loc[self.df["name"] == gene, "motifs"].values[0]

    def gene(self, id):
        try:
            return self.df.loc[id, "name"]
        except KeyError:
            return None

    def gene_to_coords(self, gene):
        return self.df.loc[self.df["name"] == gene].index[0]

    def pathogenic_min_length(self, gene):
        return self.df.loc[self.df["name"] == gene, "pathogenic_min"].values[0]

    def reflen(self, gene):
        return self.df.loc[self.df["name"] == gene, "reflen"].values[0]

    def start(self, gene):
        return self.df.loc[self.df["name"] == gene, "start"].values[0]

    def end(self, gene):
        return self.df.loc[self.df["name"] == gene, "end"].values[0]

    def chrom(self, gene):
        return self.df.loc[self.df["name"] == gene, "chrom"].values[0]

    def coords(self, gene):
        return self.df.loc[self.df["name"] == gene, ["chrom", "start", "end"]].values[0]
