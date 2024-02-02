import pandas as pd


class Repeats(object):
    def __init__(self, bed):
        self.df = self.get_repeat_info(bed)

    def get_repeat_info(self, bed):
        """
        This function parses a bed file as obtained from STRchive, which is intended for TRGT,
        but also works for STRdust.
        """
        bed = pd.read_csv(
            bed, sep="\t", header=None, names=["chrom", "start", "end", "info"]
        )
        # the TRGT bed file has as ID the <disease>_<gene> and we only care about the gene
        bed["name"] = bed["info"].apply(
            lambda x: [i.split("_")[1] for i in x.split(";") if i.startswith("ID=")][
                0
            ].replace("ID=", "")
        )
        # in case there are duplicates in the name column, use the original ID from the info field
        dups = bed.duplicated(subset="name", keep=False)
        bed.loc[dups, "name"] = bed.loc[dups, "info"].apply(
            lambda x: [
                Repeats.fix_name(i) for i in x.split(";") if i.startswith("ID=")
            ][0].replace("ID=", "")
        )
        bed["motifs"] = bed["info"].apply(
            lambda x: [
                i.replace("MOTIFS=", "").split(",")
                for i in x.split(";")
                if i.startswith("MOTIFS=")
            ][0]
        )
        bed["motif_length"] = bed["motifs"].apply(lambda x: len(x[0]))
        bed["id"] = (
            bed["chrom"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str)
        )
        return bed.drop(columns=["info", "chrom", "start", "end"]).set_index("id")

    @staticmethod
    def fix_name(name):
        return f"{name.split('_')[1]}_{name.split('_')[0]}"

    def motif_length(self, gene):
        return self.df.loc[self.df["name"] == gene, "motif_length"].values[0]

    def motifs(self, gene):
        return self.df.loc[self.df["name"] == gene, "motifs"].values[0]

    def gene(self, id):
        try:
            return self.df.loc[id, "name"]
        except KeyError:
            return None
