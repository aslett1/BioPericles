import unittest

from collections import OrderedDict
from mock import patch
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

  @patch('biopericles.TreeBuilder.tempfile.NamedTemporaryFile')
  def test_create_temporary_phylip(self, temp_mock):
    temp_file = StringIO()
    temp_mock.return_value = temp_file

    builder = TreeBuilder()

    fasta_file = StringIO("""\
>cluster_A
NNAACAAAANN
>cluster_B
CCAACAAAANN
""")

    expected_output = """\
 2 11
cluster_A  NNAACAAAAN N
cluster_B  CCAACAAAAN N
"""

    builder.load_fasta_sequences(fasta_file)
    output_file = builder._create_temporary_phylip(builder.sequences)

    self.assertEqual(output_file, temp_file)

    output_file.seek(0)
    self.assertEqual(output_file.read(), expected_output)

  def test_merge_commandline_arguments(self):
    builder = TreeBuilder()

    default_args = {'foo': 'bar'}
    new_args = {'baz': 'bar'}
    expected = {'foo': 'bar', 'baz': 'bar'}
    actual = builder._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)
    self.assertEqual(default_args, {'foo': 'bar'})

    default_args = {'foo': 'bar', 'baz': 'bar'}
    new_args = {'foo': None}
    expected = {'baz': 'bar'}
    actual = builder._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)

    default_args = {'baz': 'bar'}
    new_args = {'foo': None}
    expected = {'baz': 'bar'}
    actual = builder._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)

  def test_build_commandline_arguments(self):
    builder = TreeBuilder()

    arguments = OrderedDict([('-foo', 'bar')])
    actual = builder._build_commandline_arguments(arguments)
    self.assertEqual(actual, ['-foo', 'bar'])

    arguments = OrderedDict([('-foo', 'bar'), ('-v', ''), ('-h', None)])
    actual = builder._build_commandline_arguments(arguments)
    self.assertEqual(actual, ['-foo', 'bar', '-v'])
