from collections import Counter
import pandas as pd
import logging


def parse_kmers(df, repeats, gene):
    kmers_extracted = []
    gene_df = df[df["gene"] == gene]
    if len(gene_df) == 0:
        return pd.DataFrame()
    motif_length = repeats.motif_length(gene)
    known_motifs = repeats.motifs(gene)
    for sample, allele, seq in gene_df[["sample", "allele", "sequence"]].itertuples(
        index=False, name=None
    ):
        # check if seq is a string, as it can be a float-like NaN for reasons I don't understand
        if seq and isinstance(seq, str):
            kmers = count_kmers(seq, k=motif_length, motifs=known_motifs)
            if kmers:
                kmers.update(
                    {
                        "identifier": f"{sample}_{allele}",
                        "length": len(seq) / motif_length,
                    }
                )
                kmers_extracted.append(kmers),
    if len(kmers_extracted) == 0:
        return pd.DataFrame()
    return (
        pd.DataFrame(kmers_extracted)
        .set_index("identifier")
        .astype(float)
        .fillna(0.0)
        .round(2)
    )


def count_kmers(seq, k, motifs=None):
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmers[seq[i : i + k]] += 1
    return prune_counts(kmers, motifs)


def get_rotations(kmer, motifs):
    """
    Rotate a kmer to get all equivalent representations
    Return the lexicographical first separately, except if motifs are given
    Then prioritize known motifs. For now there is no rule on how to prioritize the known motifs.
    Also return all rotations.
    """
    e = len(kmer)
    rotations = [kmer[i:e] + kmer[:i] for i in range(e)]
    known_motifs = [m for m in motifs if m in rotations]
    if known_motifs:
        return known_motifs[0], rotations
    else:
        return sorted(rotations)[0], rotations


def prune_counts(kmers, motifs=None):
    """
    For all rotations of a kmer, keep only the lexicographical first
    Return the number as a fraction of the total kmers
    """
    pruned = dict()
    for key in kmers:
        first, rotations = get_rotations(key, motifs)
        if first in pruned.keys():
            continue
        else:
            pruned[first] = sum([kmers[r] for r in rotations])
    total_kmers = sum(pruned.values())
    return {k: v / total_kmers for k, v in pruned.items()}
