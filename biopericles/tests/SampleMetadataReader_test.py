import unittest
from mock import patch

from biopericles.SampleMetadataReader import SampleMetadataReader

class TestSampleMetadataReader(unittest.TestCase):

  @patch('biopericles.SampleMetadataReader.os')
  def test_initialise_file_not_exist(self, os_mock):
    os_mock.path.isfile.return_value = False
    self.assertRaises(ValueError, SampleMetadataReader, '/not-real')
