import unittest
from mock import patch, MagicMock
from StringIO import StringIO

from biopericles.SampleMetadataReader import SampleMetadataReader, BadMetadataException

class TestSampleMetadataReader(unittest.TestCase):

  def uninitialised_reader(self):
    # In most of the tests, we don't want to test the __init__ method
    return SampleMetadataReader.__new__(SampleMetadataReader)

  @patch('biopericles.SampleMetadataReader.os')
  def test_initialise_file_not_exist(self, os_mock):
    os_mock.path.isfile.return_value = False
    self.assertRaises(ValueError, SampleMetadataReader, '/not-real')

  def test_parse_line(self):
    reader = self.uninitialised_reader()
    reader.sample_name_idx = 1
    reader.cluster_name_idx = 0

    line = ['clusterA', 'sequence1', 'something else']
    self.assertEqual(reader.parse_line(line), ('sequence1', 'clusterA'))

    reader.sample_name_idx = 10
    self.assertRaises(IndexError, reader.parse_line, line)

  def test_create_cluster_sample_map(self):
    reader = self.uninitialised_reader()

    output = reader.create_cluster_sample_map([ ('seq_1', 'cluster_A'),
                                                ('seq_1', 'cluster_A'),
                                                ('seq_2', 'cluster_B'),
                                                ('seq_3', 'cluster_B'),
                                                ('seq_4', 'cluster_C') ])
    expected_output = {
                        'cluster_A': ['seq_1'],
                        'cluster_B': ['seq_2', 'seq_3'],
                        'cluster_C': ['seq_4']
                      }

    self.assertEqual(output, expected_output)

  def test_create_sample_cluster_map(self):
    reader = self.uninitialised_reader()

    output = reader.create_sample_cluster_map([ ('seq_1', 'cluster_A'),
                                                ('seq_1', 'cluster_A'),
                                                ('seq_2', 'cluster_B'),
                                                ('seq_3', 'cluster_B'),
                                                ('seq_4', 'cluster_C') ])
    expected_output = {
                        'seq_1': 'cluster_A',
                        'seq_2': 'cluster_B',
                        'seq_3': 'cluster_B',
                        'seq_4': 'cluster_C'
                      }

    self.assertEqual(output, expected_output)

    sample_cluster_tuples = [('seq_1', 'cluster_A'), ('seq_1', 'cluster_B')]
    self.assertRaises(BadMetadataException, reader.create_sample_cluster_map, sample_cluster_tuples)
