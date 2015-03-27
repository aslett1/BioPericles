import unittest
from mock import patch

from biopericles.ClusterSplitter import ClusterSplitter

class TestClusterSplitter(unittest.TestCase):

  def uninitialised_splitter(self):
    return ClusterSplitter.__new__(ClusterSplitter)

  def fake_abspath(self, path):
    if path[0] == '~':
      return '/home' + path[1:]
    else:
      return path

  @patch('biopericles.ClusterSplitter.os.path.abspath')
  def test_create_default_output_directory(self, abspath_mock):

    abspath_mock.side_effect = self.fake_abspath
    splitter = self.uninitialised_splitter()

    multifasta = '/parent_dir/child_dir/file.aln'
    output_directory = splitter.create_default_output_directory(multifasta)
    self.assertEqual(output_directory, '/parent_dir/child_dir')

    multifasta = '~/child_dir/file.aln'
    output_directory = splitter.create_default_output_directory(multifasta)
    self.assertEqual(output_directory, '/home/child_dir')
