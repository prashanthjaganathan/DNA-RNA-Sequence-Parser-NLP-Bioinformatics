import unittest
from dna_parser import determine_sequence_type, identify_orf, compute_sequence_stats, extract_sequence, find_type


class TestDNAParser(unittest.TestCase):
    def test_determine_sequence_type(self):
        dna_sequence = "ATGCATGCATGC"
        result = determine_sequence_type(dna_sequence)
        self.assertEqual(result, "DNA")

        # Invalid DNA
        invalid_dna_sequence = "ATGCATGCU"
        result = determine_sequence_type(invalid_dna_sequence)
        self.assertEqual(result, "Invalid")

        rna_sequence = "AUGCAUGCAUGC"
        result = determine_sequence_type(rna_sequence)
        self.assertEqual(result, "RNA")


    def test_identify_orf_dna(self):
        sequence = "ATGAAATGATT"
        orfs, count = identify_orf(sequence)
        self.assertEqual(count, 1)
        self.assertIn("ATGAAATGA", orfs)

    def test_identify_orf_rna(self):
        sequence = "AUGAAUUGA"
        orfs, count = identify_orf(sequence)
        self.assertEqual(count, 1)
        self.assertIn("AUGAAUUGA", orfs)

    def test_compute_sequence_stats_dna(self):
        sequence = "ATGAATGATAA"
        stats = compute_sequence_stats(sequence, "DNA")
        self.assertIn('Length: 11\nCounts: A: 6,C: 0,G: 2,T: 3,Ambiguous: 0,\nNumber of ORFs: 0\n', stats)

    def test_extract_sequence(self):
        # Testing the sequence extraction function
        group = ">SEQ_1\nATGAATGATAA"
        sequence = extract_sequence(group)
        self.assertEqual(sequence, "ATGAATGATAA")

    def test_find_type(self):
        # Testing the find type function
        group = ">SEQ_1\nATGAATGATAA"
        seq_type = find_type(group)
        self.assertEqual(seq_type, "DNA")

if __name__ == "__main__":
    unittest.main()