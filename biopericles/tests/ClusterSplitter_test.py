import os
import unittest
from mock import patch

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

  def fake_dirname(self, path):
    lookup = {
               '/home/another_directory': '/home',
               '/home/another_directory/': '/home/another_directory',
               '/home/child_dir': '/home',
               '/home/child_dir/': '/home/child_dir',
               '/home/child_dir/file.aln': '/home/child_dir',
               '/parent_dir/child_dir': '/parent_dir',
               '/parent_dir/child_dir/': '/parent_dir/child_dir',
               '/parent_dir/child_dir/file.aln': '/parent_dir/child_dir',
               '~/another_directory': '~',
               '~/another_directory/': '~/another_directory',
               '~/child_dir': '~',
               '~/child_dir/': '~/child_dir',
               '~/child_dir/file.aln': '~/child_dir'
             }
    return lookup[path]

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
  def test_absolute_directory_path(self, path_mock):

    path_mock.abspath.side_effect = self.fake_abspath
    path_mock.dirname.side_effect = self.fake_dirname
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

  def test_get_clusters(self):
    splitter = self.uninitialised_splitter()

    sample_to_cluster_map = {
      'seq_2': 'cluster_B',
      'seq_3': 'cluster_B',
      'seq_1': 'cluster_A',
      'seq_4': 'cluster_C'
    }

    cluster_list = splitter.get_clusters(sample_to_cluster_map)

    self.assertEqual(cluster_list, ['cluster_A', 'cluster_B', 'cluster_C'])
