import unittest
import os

from StringIO import StringIO
from mock import patch, MagicMock
from unittest import skip

from biopericles.ClusterFastaFile import ClusterFastaFile

class GeneratorMock(object):
  def __init__(self, sequences):
    self.call_count = 0
    self.sequences = (seq for seq in sequences)

  def __iter__(self):
    for seq in self.sequences:
      self.call_count += 1
      yield seq 

  def next(self):
    self.call_count += 1
    return self.sequences.next()

class TestClusterFastaFile(unittest.TestCase):
  def test_get_pair_consensus(self):
    cluster = ClusterFastaFile()

    self.assertRaises(ValueError, cluster.get_pair_consensus, 'A', 'AG')
    self.assertEqual(cluster.get_pair_consensus('AG', 'AG'), 'AG')
    self.assertEqual(cluster.get_pair_consensus('AG', 'AC'), 'AN')
    self.assertEqual(cluster.get_pair_consensus('AGT', 'ACT'), 'ANT')
    self.assertEqual(cluster.get_pair_consensus('-GT', 'ACT'), 'NNT')
    self.assertEqual(cluster.get_pair_consensus('-GT', '-CT'), 'NNT')

  def test_compare_nucleotides(self):
    cluster = ClusterFastaFile()

    self.assertEqual(cluster.compare_nucleotides('A', 'A'), 'A')
    self.assertEqual(cluster.compare_nucleotides('A', 'T'), 'N')
    self.assertEqual(cluster.compare_nucleotides('A', 'N'), 'N')
    self.assertEqual(cluster.compare_nucleotides('N', 'N'), 'N')
    self.assertEqual(cluster.compare_nucleotides('A', '-'), 'N')
    self.assertEqual(cluster.compare_nucleotides('-', '-'), 'N')

  def test_calculate_consensus(self):
    cluster = ClusterFastaFile()

    sequence_generator = (s for s in ['AG', 'AG', 'A'])
    self.assertRaises(ValueError, cluster.calculate_consensus, sequence_generator)

    sequence_generator = (s for s in ['AG', 'AG', 'AG'])
    self.assertEqual(cluster.calculate_consensus(sequence_generator), 'AG')

    sequence_generator = (s for s in ['AG', 'AG', 'TG'])
    self.assertEqual(cluster.calculate_consensus(sequence_generator), 'NG')

    sequence_generator = (s for s in ['AG', 'AT', 'TG'])
    self.assertEqual(cluster.calculate_consensus(sequence_generator), 'NN')

    sequence_generator = (s for s in ['AG', 'A-', 'AG'])
    self.assertEqual(cluster.calculate_consensus(sequence_generator), 'AN')

  @patch('biopericles.ClusterFastaFile.Bio.SeqIO')
  def test_load_file(self, seqio_mock):
    cluster = ClusterFastaFile()

    fake_file = StringIO("Some sequences")

    fake_sequence_generator = GeneratorMock(['AG', 'AG', 'AG'])

    seqio_mock.parse.return_value = fake_sequence_generator

    cluster.load_file(fake_file)
    seqio_mock.parse.assert_called_once_with(fake_file, 'fasta')
    self.assertEqual(cluster.sequences, fake_sequence_generator)
    self.assertEqual(fake_sequence_generator.call_count, 0)

    self.assertEqual(cluster.get_consensus(), 'AG')
    self.assertEqual(fake_sequence_generator.call_count, 3)

  def test_actually_load_a_file(self):
    cluster = ClusterFastaFile()
    this_dir = os.path.dirname(__file__)
    test_file_path = os.path.join(this_dir, 'data', 'cluster_A_multifasta.aln')
    with open(test_file_path, 'r') as test_file:
      cluster.load_file(test_file)
      self.assertEqual(cluster.get_consensus(), 'ANNAACAAAANN')
