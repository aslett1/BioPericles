import unittest
from mock import patch

from biopericles.SampleMetadataReader import SampleMetadataReader, BadMetadataException

class TestSampleMetadataReader(unittest.TestCase):

  @patch('biopericles.SampleMetadataReader.os')
  def test_initialise_file_not_exist(self, os_mock):
    os_mock.path.isfile.return_value = False
    self.assertRaises(ValueError, SampleMetadataReader, '/not-real')

  @patch('biopericles.SampleMetadataReader.os')
  def test_parse_line(self, os_mock):
    os_mock.path.isfile.return_value = True
    reader = SampleMetadataReader('')
    line = ['clusterA', 'sequence1', 'something else']
    self.assertEqual(reader.parse_line(line, 1, 0), ('sequence1', 'clusterA'))

    self.assertRaises(IndexError, reader.parse_line, line, 10, 0)

  @patch('biopericles.SampleMetadataReader.os')
  def test_create_cluster_sample_map(self, os_mock):
    os_mock.path.isfile.return_value = True
    reader = SampleMetadataReader('')

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

  @patch('biopericles.SampleMetadataReader.os')
  def test_create_sample_cluster_map(self, os_mock):
    os_mock.path.isfile.return_value = True
    reader = SampleMetadataReader('')

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
