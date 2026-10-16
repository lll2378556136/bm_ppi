"""Synthetic tests for Stage 1 CC-filtered fragment expansion."""

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from bm_ppi.protein_pairs import INPUT_FASTAS, run


class ProteinPairsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input_dir = self.root / "segments"
        self.input_dir.mkdir()
        self.pairs_path = self.root / "stage1.tsv.gz"
        self.output_dir = self.root / "output"
        for filename in INPUT_FASTAS:
            (self.input_dir / filename).write_text("", encoding="utf-8")

    def fasta(self, filename, text):
        (self.input_dir / filename).write_text(text, encoding="utf-8")

    def pairs(self, rows, columns=("protein_A", "protein_B", "pair_class")):
        with gzip.open(self.pairs_path, "wt", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)

    def test_filtered_cartesian_product_and_length_boundary(self):
        self.fasta("segmented.fasta", ">A-seg1\nABC\n>A-seg2\nXY\n")
        self.fasta("unsegmented.fasta", ">B-seg1\nDE\n")
        self.fasta("no_domain.fasta", ">C-seg1\nZZ\n")
        self.fasta("oversized_segments.fasta", ">D-seg1\nLONG\n")
        self.pairs(
            [
                ("A", "B", "Same locality"),
                ("A", "D", "Same locality"),
                ("A", "missing", "Same locality"),
            ]
        )
        stats = run(self.pairs_path, self.input_dir, self.output_dir, max_length=3)
        self.assertEqual(stats["fragment_pairs_written"], 2)
        self.assertEqual(stats["pairs_read"], 3)
        self.assertEqual(stats["pairs_expandable"], 1)
        self.assertEqual(stats["pairs_length_filtered"], 1)
        self.assertEqual(stats["pairs_missing_parent"], 1)
        with gzip.open(self.output_dir / "protein_pairs.fasta.gz", "rt") as handle:
            self.assertEqual(handle.read(), ">A-seg1_B-seg1\nABC:DE\n>A-seg2_B-seg1\nXY:DE\n")
        self.assertIn(
            ">A-seg1 source=segmented.fasta\nABC\n",
            (self.output_dir / "passed_sequences.fasta").read_text(),
        )
        self.assertEqual(
            (self.output_dir / "failed_sequences.fasta").read_text(),
            ">D-seg1 source=oversized_segments.fasta\nLONG\n",
        )
        with (self.output_dir / "summary.tsv").open() as handle:
            summary = dict(csv.reader(handle, delimiter="\t"))
        self.assertEqual(summary["fragment_pairs_written"], "2")
        self.assertEqual(summary["parents_failed_only"], "1")
        with self.assertRaisesRegex(ValueError, "already exist"):
            run(self.pairs_path, self.input_dir, self.output_dir, 3)

    def test_all_invalid_inputs_rejected(self):
        self.fasta("segmented.fasta", ">A-seg1\nAC\n")
        cases = [
            (">A-seg1\nAC\n>A-seg1\nGG\n", "Duplicate segment ID"),
            (">A\nAC\n", "Invalid segment ID"),
            (">A-seg1\n", "Empty or invalid sequence"),
            (">A-seg0\nAC\n", "Invalid segment ID"),
            ("not a fasta\n", "Invalid FASTA header"),
        ]
        self.pairs([("A", "B", "Same locality")])
        for content, error in cases:
            with self.subTest(error=error):
                self.fasta("segmented.fasta", content)
                with self.assertRaisesRegex(ValueError, error):
                    run(self.pairs_path, self.input_dir, self.output_dir)
                self.assertFalse(self.output_dir.exists())

        self.fasta("segmented.fasta", ">A-seg1\nAC\n")
        self.pairs([("A", "B", "Same locality")], columns=("protein_A", "protein_B"))
        with self.assertRaisesRegex(ValueError, "columns"):
            run(self.pairs_path, self.input_dir, self.output_dir)
        self.pairs([("A", "A", "Same locality")])
        with self.assertRaisesRegex(ValueError, "Invalid Stage 1 pair"):
            run(self.pairs_path, self.input_dir, self.output_dir)

    def test_bad_pair_class_and_row_width(self):
        self.fasta("segmented.fasta", ">A-seg1\nAA\n")
        for rows in [
            [("A", "B", "Different")],
            [("A", "B")],
            [("A", "B", "Same locality", "extra")],
        ]:
            with self.subTest(rows=rows):
                self.pairs(rows)
                with self.assertRaisesRegex(ValueError, "Stage 1"):
                    run(self.pairs_path, self.input_dir, self.root / str(len(rows[0])))

    def test_missing_fasta_and_nonpositive_length(self):
        self.pairs([])
        with self.assertRaisesRegex(ValueError, "max_length"):
            run(self.pairs_path, self.input_dir, self.output_dir, 0)
        (self.input_dir / "no_domain.fasta").unlink()
        with self.assertRaisesRegex(ValueError, "Missing required"):
            run(self.pairs_path, self.input_dir, self.output_dir)


if __name__ == "__main__":
    unittest.main()
