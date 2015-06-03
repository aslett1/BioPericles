import unittest

from collections import OrderedDict
from StringIO import StringIO

from biopericles.Common import LoadFastaMixin, RunExternalApplicationMixin

class TestLoadFastaMixin(unittest.TestCase):
  def test_load_fasta_sequences(self):
    mixin = LoadFastaMixin()

    fasta_file = StringIO("""\
>cluster_A
NNAACAAAANN
>cluster_B
CCAACAAAANN
""")

    mixin.load_fasta_sequences(fasta_file)

    self.assertItemsEqual(mixin.sequences.keys(), ['cluster_A', 'cluster_B'])
    self.assertItemsEqual(mixin.sequences['cluster_A'].seq, "NNAACAAAANN")
    self.assertItemsEqual(mixin.sequences['cluster_B'].seq, "CCAACAAAANN")

    fasta_file = StringIO("""\
>cluster_A
GGAACAAAANN
""")

    mixin.load_fasta_sequences(fasta_file)

    self.assertItemsEqual(mixin.sequences.keys(), ['cluster_A'])
    self.assertItemsEqual(mixin.sequences['cluster_A'].seq, "GGAACAAAANN")

    fasta_file = StringIO("""\
>cluster_X
GGAACAAAANN
>cluster_X
TTAACAAAANN
""")

    self.assertRaises(ValueError, mixin.load_fasta_sequences, fasta_file)

class TestRunExternalApplicationMixin(unittest.TestCase):
  def test_run_application(self):
    # Tested pretty well by TreeBuilder
    pass

  def test_merge_commandline_arguments(self):
    mixin = RunExternalApplicationMixin()

    default_args = {'foo': 'bar'}
    new_args = {'baz': 'bar'}
    expected = {'foo': 'bar', 'baz': 'bar'}
    actual = mixin._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)
    self.assertEqual(default_args, {'foo': 'bar'})

    default_args = {'foo': 'bar', 'baz': 'bar'}
    new_args = {'foo': None}
    expected = {'baz': 'bar'}
    actual = mixin._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)

    default_args = {'baz': 'bar'}
    new_args = {'foo': None}
    expected = {'baz': 'bar'}
    actual = mixin._merge_commandline_arguments(default_args, new_args)
    self.assertEqual(expected, actual)

  def test_build_commandline_arguments(self):
    mixin = RunExternalApplicationMixin()

    arguments = OrderedDict([('-foo', 'bar')])
    actual = mixin._build_commandline_arguments(arguments)
    self.assertEqual(actual, ['-foo', 'bar'])

    arguments = OrderedDict([('-foo', 'bar'), ('-v', ''), ('-h', None)])
    actual = mixin._build_commandline_arguments(arguments)
    self.assertEqual(actual, ['-foo', 'bar', '-v'])
