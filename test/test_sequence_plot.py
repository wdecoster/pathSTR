import unittest
from types import SimpleNamespace

import pandas as pd

from pathSTR.plot import kmer_plot_sequence
from scripts.aSTRonaut import plot_sequence


class SequencePlotTests(unittest.TestCase):
    def setUp(self):
        self.repeats = pd.DataFrame(
            {
                "sample": ["short", "long", "medium"],
                "allele": ["1", "2", "1"],
                "sequence": ["AT", "ATATAT", "ATAT"],
                "length": [2, 6, 4],
                "case": [0, 0, 0],
            },
            index=[10, 20, 30],
        )
        self.kmers = pd.DataFrame({"AT": [1, 3, 2], "length": [2, 6, 4]})

    def test_app_order_with_filtered_indices_and_direction(self):
        for direction in ("left-to-right", "right-to-left"):
            with self.subTest(direction=direction):
                fig = kmer_plot_sequence(
                    self.kmers,
                    self.repeats,
                    "test",
                    length_range=[3, 6],
                    direction=direction,
                )
                self.assertEqual(
                    tuple(fig.layout.yaxis.categoryarray), ("long_2", "medium_1")
                )
                self.assertEqual(
                    sum(len(trace.x) for trace in fig.data), 10
                )
                if direction == "right-to-left":
                    self.assertEqual(fig.layout.xaxis.autorange, "reversed")

    def test_astronaut_order_with_non_default_indices(self):
        args = SimpleNamespace(
            alphabetic=False,
            hide_allele_label=False,
            somatic=False,
            size=3,
            hide_labels=False,
            label_size=10,
            publication=False,
            minimal=False,
            height=800,
            width=None,
            legend_corner="outside",
        )
        fig = plot_sequence(self.repeats.copy(), ["AT"], "chr1:1", args)
        self.assertEqual(
            tuple(fig.layout.yaxis.categoryarray), ("short_1", "medium_1", "long_2")
        )
        self.assertEqual(sum(len(trace.x) for trace in fig.data), 12)


if __name__ == "__main__":
    unittest.main()
