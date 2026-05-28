import pandas as pd
import logging


class Repeats(object):
    def __init__(self, df=None, build=["hg38", "t2t"]):
        if df is not None:
            self.df = df
        else:
            self.df = pd.concat(self.get_repeat_info(b) for b in build)
        self._build_overlap_index()

    def _build_overlap_index(self):
        """Per-(build, chrom) list of (start, end, name) for overlap-based locus lookup.
        Used by `coords_to_gene` because the (start, end) reported by callers can vary by a
        few bp depending on the genotyper's anchoring convention (STRdust historically used
        POS=BED.start, then POS=BED.start+1, current builds POS=BED.start+2). An exact
        coordinate match is too brittle; an overlap test tolerates all conventions while
        STRchive's loci are far enough apart that ambiguity is not a concern."""
        self._intervals = {}
        for idx, row in self.df.iterrows():
            chrom, coords = idx.split(":")
            start_s, end_s = coords.split("-")
            self._intervals.setdefault((row["build"], chrom), []).append(
                (int(start_s), int(end_s), row["name"])
            )

    def get_repeat_info(self, build):
        """
        This function parses a bed file as obtained from STRchive, which is intended for TRGT, but also works for STRdust.
        Also a CSV is downloaded, with extra information about the repeats.
        A dataframe is constructed, which is the data underlying the Repeats class instance

        The motif length that is used for e.g. kmer plots is the shortest motif that is not 1.
        Using the longest motif causes kmer file sizes to explode for genes with long motifs.
        """
        urls = {
            "hg38": "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.TRGT.bed",
            "t2t": "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.T2T-chm13.TRGT.bed",
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
            pd.read_json(
                "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json").set_index("id")[["pathogenic_min"]],
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
        bed["motif_length"] = bed["motifs"].apply(self.get_motif_length)
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
    def get_motif_length(motifs):
        """
        Return the length of the shortest motif that is not 1
        Falls back to 1 if all motifs have length 1
        """
        lengths = sorted([len(m) for m in motifs])
        # Return the first motif length that is not 1, or 1 if all are 1
        for length in lengths:
            if length > 1:
                return length
        return lengths[0]

    def strchive_version(self):
        """get the latest release from https://github.com/dashnowlab/STRchive"""
        import requests

        url = "https://api.github.com/repos/dashnowlab/STRchive/releases/latest"
        response = requests.get(url)
        return response.json()["tag_name"]

    def query(self, dataset, gene, column):
        """
        Generic method to query the dataframe based on a dataset, gene and column (or list of columns).
        Raises KeyError if no row matches — the caller is responsible for handling it. We don't
        sys.exit here because this method is called from Dash callbacks; a missing gene must not
        bring down the whole server.
        """
        build = dataset.split("_")[1]
        try:
            return self.df.loc[
                (self.df["build"] == build) & (self.df["name"] == gene), column
            ].values[0]
        except IndexError as e:
            logging.error(f"Query in repeats failed for build={build}, gene={gene}, column={column}")
            raise KeyError(
                f"No repeat entry for build={build}, gene={gene}"
            ) from e

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

    def coords_to_gene(self, chrom, start, end, build):
        """Return the catalog gene whose interval overlaps the query, or None.

        Overlap is half-open: a catalog [s, e) and query [start, end) match iff s < end
        and e > start. This tolerates the small anchor/padding differences that arise
        between STRdust versions (and between STRdust and LongTR) — e.g. POS = BED.start
        + 0, +1, or +2 all still resolve to the same locus.
        """
        chrom_norm = chrom if chrom.startswith("chr") else f"chr{chrom}"
        intervals = self._intervals.get((build, chrom_norm), [])
        matches = [name for (s, e, name) in intervals if s < end and e > start]
        if not matches:
            return None
        if len(matches) > 1:
            logging.warning(
                f"Multiple catalog hits for {chrom_norm}:{start}-{end} in {build}: "
                f"{matches} — returning first."
            )
        return matches[0]

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
