from collections import Counter
import pandas as pd
import logging


def parse_kmers(df, repeats, gene):
    kmers_extracted = []
    gene_df = df[df["gene"] == gene]
    motif_length = repeats.motif_length(gene)
    for sample, allele, seq in gene_df[["sample", "allele", "sequence"]].itertuples(
        index=False, name=None
    ):
        if seq:
            kmers = count_kmers(seq, k=motif_length)
            if kmers:
                kmers.update({"identifier": f"{sample}_{allele}", "length": len(seq)})
                kmers_extracted.append(kmers),
    return (
        pd.DataFrame(kmers_extracted)
        .set_index("identifier")
        .astype(float)
        .fillna(0.0)
        .round(2)
    )


def count_kmers(seq, k):
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmers[seq[i : i + k]] += 1
    return prune_counts(kmers)


def get_rotations(kmer):
    """
    Rotate a kmer to get all equivalent representations
    """
    e = len(kmer)
    rotations = [kmer[i:e] + kmer[:i] for i in range(e)]
    return sorted(rotations)[0], rotations


def prune_counts(kmers):
    """
    For all rotations of a kmer, keep only the lexicographical first
    Return the number as a fraction of the total kmers
    """
    pruned = dict()
    for key in kmers:
        first, rotations = get_rotations(key)
        if first in pruned.keys():
            continue
        else:
            pruned[first] = sum([kmers[r] for r in rotations])
    total_kmers = sum(pruned.values())
    return {k: v / total_kmers for k, v in pruned.items()}
