import os
import unittest
from mock import patch
from StringIO import StringIO

from biopericles.ClusterSplitter import ClusterSplitter
from biopericles.ClusterSplitter import NotDirectoryException

class TestClusterSplitter(unittest.TestCase):

  def uninitialised_splitter(self):
    return ClusterSplitter.__new__(ClusterSplitter)

  def fake_abspath(self, path):
    if path[0] == '~':
      return '/home' + path[1:]
    else:
      return path

  def fake_normpath(self, path):
    lookup = {
               '/home/another_directory': '/home/another_directory',
               '/home/another_directory/': '/home/another_directory',
               '/home/child_dir': '/home/child_dir',
               '/home/child_dir/': '/home/child_dir',
               '/home/child_dir/file.aln': '/home/child_dir/file.aln',
               '/parent_dir/child_dir': '/parent_dir/child_dir',
               '/parent_dir/child_dir/': '/parent_dir/child_dir',
               '/parent_dir/child_dir/file.aln': '/parent_dir/child_dir/file.aln',
               '~/another_directory': '~/another_directory',
               '~/another_directory/': '~/another_directory',
               '~/child_dir': '~/child_dir',
               '~/child_dir/': '~/child_dir',
               '~/child_dir/file.aln': '~/child_dir/file.aln'
             }
    return lookup[path]

  @patch('biopericles.ClusterSplitter.ClusterSplitter.absolute_directory_path')
  @patch('biopericles.ClusterSplitter.os.getcwd')
  def test_init_output_directory(self, cwd_mock, directory_path_mock):
    cwd_mock.return_value = '/parent_dir/child_dir'

    multifasta = '/parent_dir/child_dir/file.aln'
    sequence_to_cluster_map = {}

    directory_path_mock.return_value = '/parent_dir/child_dir'
    splitter = ClusterSplitter(multifasta, sequence_to_cluster_map)
    self.assertEqual(splitter.output_directory, '/parent_dir/child_dir')

    directory_path_mock.return_value = '/home/another_directory'
    splitter = ClusterSplitter(multifasta, sequence_to_cluster_map, '~/another_directory/')
    self.assertEqual(splitter.output_directory, '/home/another_directory')


  @patch('biopericles.ClusterSplitter.os.path')
  def test_absolute_multifasta_path(self, path_mock):

    path_mock.abspath.side_effect = self.fake_abspath
    path_mock.normpath.side_effect = self.fake_normpath
    path_mock.sep.return_value = '/'
    splitter = self.uninitialised_splitter()

    path_mock.isdir.return_value = False

    multifasta_input_path = '/parent_dir/child_dir/file.aln'
    self.assertEqual(splitter.absolute_multifasta_path(multifasta_input_path), '/parent_dir/child_dir/file.aln')

    multifasta_input_path = '~/child_dir/file.aln'
    self.assertEqual(splitter.absolute_multifasta_path(multifasta_input_path), '/home/child_dir/file.aln')

    path_mock.isdir.return_value = True

    directory = '~/another_directory/'
    self.assertRaises(NotDirectoryException, splitter.absolute_directory_path, )

    directory = '~/another_directory'
    output_directory = splitter.absolute_directory_path(directory)
    self.assertEqual(output_directory, '/home/another_directory')

  @patch('biopericles.ClusterSplitter.os.path')
  def test_absolute_directory_path(self, path_mock):

    path_mock.abspath.side_effect = self.fake_abspath
    path_mock.normpath.side_effect = self.fake_normpath
    path_mock.sep.return_value = '/'
    splitter = self.uninitialised_splitter()

    path_mock.isdir.return_value = False

    multifasta = '/parent_dir/child_dir/file.aln'
    self.assertRaises(NotDirectoryException, splitter.absolute_directory_path, multifasta)

    multifasta = '~/child_dir/file.aln'
    self.assertRaises(NotDirectoryException, splitter.absolute_directory_path, multifasta)

    path_mock.isdir.return_value = True

    directory = '~/another_directory/'
    output_directory = splitter.absolute_directory_path(directory)
    self.assertEqual(output_directory, '/home/another_directory')

    directory = '~/another_directory'
    output_directory = splitter.absolute_directory_path(directory)
    self.assertEqual(output_directory, '/home/another_directory')

  @patch('biopericles.ClusterSplitter.open', create=True)
  def test_create_cluster_output_files(self, open_mock):

    splitter = self.uninitialised_splitter()
    splitter.output_directory = '/home/output'

    a_file = StringIO('A file')
    open_mock.return_value = a_file

    clusters = ['cluster_A', 'cluster_B', 'cluster_C']
    splitter.create_cluster_output_files(clusters)

    self.assertItemsEqual(splitter.cluster_output_files.keys(), ['cluster_A', 'cluster_B', 'cluster_C'])
    for cluster, file_handler in splitter.cluster_output_files.items():
      self.assertEqual(file_handler, a_file, "File handler for %s isn't a file handler" % cluster)

    open_mock.assert_any_call('/home/output/cluster_cluster_A_multifasta.aln', 'w')
    open_mock.assert_any_call('/home/output/cluster_cluster_B_multifasta.aln', 'w')
    open_mock.assert_any_call('/home/output/cluster_cluster_C_multifasta.aln', 'w')

    splitter.cluster_output_files = {}

    open_mock.side_effect = IOError("Problem")
    self.assertRaises(IOError, splitter.create_cluster_output_files, clusters)

  def test_get_clusters(self):
    splitter = self.uninitialised_splitter()

    sample_to_cluster_map = {
      'seq2': 'cluster_B',
      'seq3': 'cluster_B',
      'seq1': 'cluster_A',
      'seq4': 'cluster_C'
    }

    cluster_list = splitter.get_clusters(sample_to_cluster_map)

    self.assertEqual(cluster_list, ['cluster_A', 'cluster_B', 'cluster_C'])

  def test_get_sequences(self):
    splitter = self.uninitialised_splitter()

    sample_to_cluster_map = {
      'seq2': 'cluster_B',
      'seq3': 'cluster_B',
      'seq1': 'cluster_A',
      'seq4': 'cluster_C'
    }

    sequence_list = splitter.get_sequences(sample_to_cluster_map)

    self.assertEqual(sequence_list, ['seq1', 'seq2', 'seq3', 'seq4'])
