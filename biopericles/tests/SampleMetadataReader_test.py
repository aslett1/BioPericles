import unittest
from mock import patch

from biopericles.SampleMetadataReader import SampleMetadataReader

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
