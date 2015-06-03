import os
import unittest

from StringIO import StringIO
from mock import patch, MagicMock
from unittest import skip

from biopericles.ClusterConsensus import ClusterConsensus

class GeneratorMock(object):
  def __init__(self, sequences):
    self.call_count = 0
    self.sequences = (FakeSequence(seq) for seq in sequences)

  def __iter__(self):
    for seq in self.sequences:
      self.call_count += 1
      yield seq 

  def next(self):
    self.call_count += 1
    return self.sequences.next()

class FakeSequence(object):
  def __init__(self, sequence):
    self.seq = sequence

class TestClusterConsensus(unittest.TestCase):
  def test_calculate_consensus(self):
    cluster = ClusterConsensus()

    sequence_generator = (FakeSequence(s) for s in ['AG', 'AG', 'A'])
    self.assertRaises(ValueError, cluster._calculate_consensus, sequence_generator)

    sequence_generator = (FakeSequence(s) for s in ['AG', 'AG', 'AG'])
    self.assertEqual(cluster._calculate_consensus(sequence_generator), 'AG')

    sequence_generator = (FakeSequence(s) for s in ['AG', 'AG', 'TG'])
    self.assertEqual(cluster._calculate_consensus(sequence_generator), 'NG')

    sequence_generator = (FakeSequence(s) for s in ['AG', 'AT', 'TG'])
    self.assertEqual(cluster._calculate_consensus(sequence_generator), 'NN')

    sequence_generator = (FakeSequence(s) for s in ['AG', 'A-', 'AG'])
    self.assertEqual(cluster._calculate_consensus(sequence_generator), 'AN')

  @patch('biopericles.ClusterConsensus.Bio.SeqIO')
  def test_load_fasta_file(self, seqio_mock):
    cluster = ClusterConsensus()

    fake_file = StringIO("Some sequences")

    fake_sequence_generator = GeneratorMock(['AG', 'AG', 'AG'])

    seqio_mock.parse.return_value = fake_sequence_generator

    cluster.load_fasta_file(fake_file)
    seqio_mock.parse.assert_called_once_with(fake_file, 'fasta')
    self.assertEqual(cluster.sequences, fake_sequence_generator)
    self.assertEqual(fake_sequence_generator.call_count, 0)

    self.assertEqual(cluster.get_consensus(), 'AG')
    self.assertEqual(fake_sequence_generator.call_count, 3)

  @patch('biopericles.ClusterConsensus.random.randrange')
  def test_get_cluster_name(self, random_mock):
    cluster = ClusterConsensus()

    cluster.output_file = StringIO()
    random_mock.return_value = 10000
    self.assertEqual(cluster._get_cluster_name(), 'cluster_10000')
    
    cluster.output_file = MagicMock()
    cluster.output_file.name = 'foo/bar.baz.mfa'
    self.assertEqual(cluster._get_cluster_name(), 'bar.baz')

    cluster.cluster_name = 'foo'
    self.assertEqual(cluster._get_cluster_name(), 'foo')

  def test_write_consensus(self):
    cluster = ClusterConsensus()

    cluster.output_file = StringIO()
    cluster.cluster_name = 'cluster_1000'
    cluster._calculate_consensus = MagicMock(return_value='NNANTGUNCN')

    cluster.write_consensus()

    expected_results = """\
>cluster_1000
NNANTGUNCN
"""
    cluster.output_file.seek(0)
    self.assertEqual(cluster.output_file.read(), expected_results)

  def test_real_file(self):
    test_folder = os.path.abspath(os.path.dirname(__file__))
    test_data_path = os.path.join(test_folder, 'data', 'cluster_A_multifasta.aln')

    cluster = ClusterConsensus()
    cluster.output_file = StringIO()
    cluster.output_file.name = "fake_filename.aln"
    with open(test_data_path, 'r') as input_file:
      cluster.load_fasta_file(input_file)
      cluster.write_consensus()

    expected_response = """\
>fake_filename
ANNAACAAAANN
"""
    cluster.output_file.seek(0)
    self.assertEqual(cluster.output_file.read(), expected_response)

  def test_get_all_sequences(self):
    test_folder = os.path.abspath(os.path.dirname(__file__))
    test_data_path = os.path.join(test_folder, 'data', 'cluster_A_multifasta.aln')

    cluster = ClusterConsensus()
    with open(test_data_path, 'r') as input_file:
      cluster.load_fasta_file(input_file)
      sequences = cluster._get_all_sequences(cluster.sequences)

    expected = ["AAAAACAAAAA-", "ATGAACAAAA--"]
    self.assertEqual(sequences, expected)

  def test_get_length_of_longest(self):
    cluster = ClusterConsensus()
    sequence_strings = ["AAAAA", "ATGAACAAAA--"]
    longest = cluster._get_length_of_longest(sequence_strings)

    self.assertEqual(longest, 12)

  def test_compare_bases(self):
    cluster = ClusterConsensus()

    self.assertEqual(cluster._compare_bases(['A', 'A']), 'A')
    self.assertEqual(cluster._compare_bases(['A', 'T']), 'N')
    self.assertEqual(cluster._compare_bases(['A', 'N']), 'N')
    self.assertEqual(cluster._compare_bases(['N', 'N']), 'N')
    self.assertEqual(cluster._compare_bases(['A', '-']), 'N')
    self.assertEqual(cluster._compare_bases(['-', '-']), 'N')
