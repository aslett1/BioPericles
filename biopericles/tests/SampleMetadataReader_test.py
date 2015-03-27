import os
import unittest
from mock import patch, MagicMock
from StringIO import StringIO

from biopericles.SampleMetadataReader import SampleMetadataReader, BadMetadataException

class TestSampleMetadataReader(unittest.TestCase):

  def uninitialised_reader(self):
    # In most of the tests, we don't want to test the __init__ method
    return SampleMetadataReader.__new__(SampleMetadataReader)

  def test_initialise_from_file(self):
    test_dir = os.path.dirname(os.path.abspath(__file__))
    path_to_csv = os.path.join(test_dir, 'data', 'clusters_spreadsheet.csv')
    reader = SampleMetadataReader(path_to_csv, 0, 2)

    expected_cluster_sample_map = {
                                     '1': ['seq1', 'seq2'],
                                     '2': ['seq3', 'seq4', 'seq5'],
                                     '3': ['seq6', 'seq7', 'seq8'],
                                     '4': ['seq10', 'seq9']
                                  }

    self.assertEqual(reader.cluster_sample_map, expected_cluster_sample_map)

    expected_sample_cluster_map = {
                                    'seq1': '1',
                                    'seq2': '1',
                                    'seq3': '2',
                                    'seq4': '2',
                                    'seq5': '2',
                                    'seq6': '3',
                                    'seq7': '3',
                                    'seq8': '3',
                                    'seq9': '4',
                                    'seq10': '4'
                                  }

    self.assertEqual(reader.sample_cluster_map, expected_sample_cluster_map)

  @patch('biopericles.SampleMetadataReader.os')
  def test_initialise_file_not_exist(self, os_mock):
    os_mock.path.isfile.return_value = False
    not_used = 0
    self.assertRaises(ValueError, SampleMetadataReader, '/not-real', not_used, not_used)

  @patch('biopericles.SampleMetadataReader.SampleMetadataReader.config_file_exists')
  @patch('biopericles.SampleMetadataReader.open', create=True)
  def test_initialise(self, open_mock, config_mock):
    fake_file_contents = """\
Sample,Cluster,Year
seq_1,cluster_A,2011
seq_1,cluster_A,2012
seq_2,cluster_B,2013
seq_3,cluster_B,2014
seq_4,cluster_C,2010
"""
    fake_file = StringIO(fake_file_contents)

    open_mock.return_value.__enter__.return_value = fake_file
    config_mock.return_value = True

    reader = SampleMetadataReader(metadata_filename='/some-file', sample_name_idx=0, cluster_name_idx=1)

    expected_cluster_sample_map = {
                        'cluster_A': ['seq_1'],
                        'cluster_B': ['seq_2', 'seq_3'],
                        'cluster_C': ['seq_4']
                      }

    self.assertEqual(reader.cluster_sample_map, expected_cluster_sample_map)

    expected_sample_cluster_map = {
                        'seq_1': 'cluster_A',
                        'seq_2': 'cluster_B',
                        'seq_3': 'cluster_B',
                        'seq_4': 'cluster_C'
                      }

    self.assertEqual(reader.sample_cluster_map, expected_sample_cluster_map)

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

  def test_parse_metadata_file(self):
    reader = self.uninitialised_reader()
    reader.sample_name_idx = 0
    reader.cluster_name_idx = 1

    fake_file_contents = """\
Sample,Cluster,Year
seq_1,cluster_A,2011
seq_1,cluster_A,2012
seq_2,cluster_B,2013
seq_3,cluster_B,2014
seq_4,cluster_C,2010
"""
    fake_file = StringIO(fake_file_contents)

    reader.parse(fake_file)

    expected_cluster_sample_map = {
                        'cluster_A': ['seq_1'],
                        'cluster_B': ['seq_2', 'seq_3'],
                        'cluster_C': ['seq_4']
                      }

    self.assertEqual(reader.cluster_sample_map, expected_cluster_sample_map)

    expected_sample_cluster_map = {
                        'seq_1': 'cluster_A',
                        'seq_2': 'cluster_B',
                        'seq_3': 'cluster_B',
                        'seq_4': 'cluster_C'
                      }

    self.assertEqual(reader.sample_cluster_map, expected_sample_cluster_map)
