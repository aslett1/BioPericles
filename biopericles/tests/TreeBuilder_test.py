import unittest

from StringIO import StringIO

from biopericles.TreeBuilder import TreeBuilder

class TestTreeBuilder(unittest.TestCase):
  def test_load_fasta_sequences(self):
    builder = TreeBuilder()

    fasta_file = StringIO("""\
>cluster_A
NNAACAAAANN
>cluster_B
CCAACAAAANN
""")

    builder.load_fasta_sequences(fasta_file)

    self.assertItemsEqual(builder.sequences.keys(), ['cluster_A', 'cluster_B'])
    self.assertItemsEqual(builder.sequences['cluster_A'].seq, "NNAACAAAANN")
    self.assertItemsEqual(builder.sequences['cluster_B'].seq, "CCAACAAAANN")

    fasta_file = StringIO("""\
>cluster_A
GGAACAAAANN
""")

    builder.load_fasta_sequences(fasta_file)

    self.assertItemsEqual(builder.sequences.keys(), ['cluster_A'])
    self.assertItemsEqual(builder.sequences['cluster_A'].seq, "GGAACAAAANN")

    fasta_file = StringIO("""\
>cluster_X
GGAACAAAANN
>cluster_X
TTAACAAAANN
""")

    self.assertRaises(ValueError, builder.load_fasta_sequences, fasta_file)
