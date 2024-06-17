from itertools import groupby


def rle(sequence, motif_length):
    """
    This function takes a sequence and a motif length and returns a compact run-length encoded sequence
    Note that this is currently not the shortest possible solution, but it is a good starting point
    Shorter could be achieved by skipping nucleotides that are not part of the motif, i.e. leaving insertions
    This would require much more complex hocus pocus (dag dynamic programming) that would eventually be implemented into STRdust

    kudos to https://stackoverflow.com/a/78634539/6631639
    """
    solutions = []
    for start in range(0, motif_length):
        parts = [
            sequence[i : i + motif_length]
            for i in range(start, len(sequence), motif_length)
        ]
        lst = [sequence[0:start] if start > 0 else ""]
        for k, g in groupby(parts):
            _len = len(list(g))
            if _len > 1:
                lst.append(f"({k}){_len}")
            else:
                lst.append(k)
        solutions.append("".join(lst))
    # return the shortest solution
    return min(solutions, key=len)
