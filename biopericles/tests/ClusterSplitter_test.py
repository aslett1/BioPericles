import os
import unittest
from mock import patch, MagicMock
from StringIO import StringIO

from biopericles.ClusterSplitter import ClusterSplitter
from biopericles.ClusterSplitter import NotFileException, NotDirectoryException, \
                                        OutputFileAlreadyExistsException

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

  def fake_sequence(self, number, base):
    seq = MagicMock()
    name = "seq%s" % number
    seq.id = name
    seq.format.return_value = ">{name}\n{sequence}\n".format(name=name, sequence=base*8)
    return seq

  @patch('biopericles.ClusterSplitter.ClusterSplitter.absolute_directory_path')
  @patch('biopericles.ClusterSplitter.ClusterSplitter.absolute_multifasta_path')
  @patch('biopericles.ClusterSplitter.os.getcwd')
  def test_init_multifasta_path(self, cwd_mock, multifasta_path_mock, directory_path_mock):
    cwd_mock.return_value = '/parent_dir/child_dir'

    sequence_cluster_map = {}

    multifasta = '/parent_dir/child_dir/file.aln'
    multifasta_path_mock.return_value = '/parent_dir/child_dir/file.aln'
    splitter = ClusterSplitter(multifasta, sequence_cluster_map)
    self.assertEqual(splitter.multifasta_path, '/parent_dir/child_dir/file.aln')

    multifasta = '~/another_directory/file.aln'
    multifasta_path_mock.return_value = '/home/another_directory/file.aln'
    splitter = ClusterSplitter(multifasta, sequence_cluster_map)
    self.assertEqual(splitter.multifasta_path, '/home/another_directory/file.aln')

  @patch('biopericles.ClusterSplitter.ClusterSplitter.absolute_directory_path')
  @patch('biopericles.ClusterSplitter.ClusterSplitter.absolute_multifasta_path')
  @patch('biopericles.ClusterSplitter.os.getcwd')
  def test_init_output_directory(self, cwd_mock, multifasta_path_mock, directory_path_mock):
    cwd_mock.return_value = '/parent_dir/child_dir'

    multifasta = '/parent_dir/child_dir/file.aln'
    sequence_cluster_map = {}

    directory_path_mock.return_value = '/parent_dir/child_dir'
    splitter = ClusterSplitter(multifasta, sequence_cluster_map)
    self.assertEqual(splitter.output_directory, '/parent_dir/child_dir')

    directory_path_mock.return_value = '/home/another_directory'
    splitter = ClusterSplitter(multifasta, sequence_cluster_map, '~/another_directory/')
    self.assertEqual(splitter.output_directory, '/home/another_directory')

  @patch('biopericles.ClusterSplitter.os.path')
  def test_absolute_multifasta_path(self, path_mock):

    path_mock.abspath.side_effect = self.fake_abspath
    path_mock.normpath.side_effect = self.fake_normpath
    path_mock.sep.return_value = '/'
    splitter = self.uninitialised_splitter()

    path_mock.isfile.return_value = True

    multifasta_input_path = '/parent_dir/child_dir/file.aln'
    self.assertEqual(splitter.absolute_multifasta_path(multifasta_input_path), '/parent_dir/child_dir/file.aln')

    multifasta_input_path = '~/child_dir/file.aln'
    self.assertEqual(splitter.absolute_multifasta_path(multifasta_input_path), '/home/child_dir/file.aln')

    path_mock.isfile.return_value = False

    directory = '~/another_directory/'
    self.assertRaises(NotFileException, splitter.absolute_multifasta_path, directory)

    directory = '~/another_directory'
    self.assertRaises(NotFileException, splitter.absolute_multifasta_path, directory)

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

  @patch('biopericles.ClusterSplitter.os.path')
  @patch('biopericles.ClusterSplitter.open', create=True)
  def test_create_cluster_output_files(self, open_mock, path_mock):

    splitter = self.uninitialised_splitter()
    splitter.output_directory = '/home/output'

    a_file = StringIO('A file')
    open_mock.return_value = a_file

    path_mock.isfile.return_value = False
    path_mock.join.side_effect = lambda *paths: "/".join(paths)

    clusters = ['cluster_A', 'cluster_B', 'cluster_C']
    splitter.create_cluster_output_files(clusters)

    self.assertItemsEqual(splitter.cluster_output_files.keys(), ['cluster_A', 'cluster_B', 'cluster_C'])
    for cluster, file_handler in splitter.cluster_output_files.items():
      self.assertEqual(file_handler, a_file, "File handler for %s isn't a file handler" % cluster)

    open_mock.assert_any_call('/home/output/cluster_cluster_A.aln', 'w')
    open_mock.assert_any_call('/home/output/cluster_cluster_B.aln', 'w')
    open_mock.assert_any_call('/home/output/cluster_cluster_C.aln', 'w')

    splitter.cluster_output_files = {}

    open_mock.side_effect = IOError("Problem")
    self.assertRaises(IOError, splitter.create_cluster_output_files, clusters)

  @patch('biopericles.ClusterSplitter.os.path')
  @patch('biopericles.ClusterSplitter.open', create=True)
  def test_create_cluster_output_files_already_exists(self, open_mock, path_mock):

    splitter = self.uninitialised_splitter()
    splitter.output_directory = '/home/output'

    path_mock.isfile.return_value = True
    path_mock.join.side_effect = lambda *paths: "/".join(paths)

    clusters = ['cluster_A']
    self.assertRaises(OutputFileAlreadyExistsException,
                      splitter.create_cluster_output_files, clusters)

    self.assertEqual(open_mock.call_count, 0)

  def test_get_clusters(self):
    splitter = self.uninitialised_splitter()

    sample_cluster_map = {
      'seq2': 'cluster_B',
      'seq3': 'cluster_B',
      'seq1': 'cluster_A',
      'seq4': 'cluster_C'
    }

    cluster_list = splitter.get_clusters(sample_cluster_map)

    self.assertEqual(cluster_list, ['cluster_A', 'cluster_B', 'cluster_C'])

  def test_get_sequences(self):
    splitter = self.uninitialised_splitter()

    sample_cluster_map = {
      'seq2': 'cluster_B',
      'seq3': 'cluster_B',
      'seq1': 'cluster_A',
      'seq4': 'cluster_C'
    }

    sequence_list = splitter.get_sequences(sample_cluster_map)

    self.assertEqual(sequence_list, ['seq1', 'seq2', 'seq3', 'seq4'])

  def test_write_sequence_to_cluster(self):
    seq = MagicMock()

    splitter = self.uninitialised_splitter()
    splitter.cluster_output_files = {'cluster_A': StringIO(), 'cluster_B': StringIO()}

    sequence_cluster_map = {'seq1': 'cluster_A', 'seq2': 'cluster_A', 'seq3': 'cluster_B'}

    seq = self.fake_sequence(1, 'A')
    success = splitter.write_sequence_to_cluster(sequence_cluster_map, seq)
    self.assertTrue(success)
    result_file = splitter.cluster_output_files['cluster_A']
    result_file.seek(0)
    result = result_file.read()
    self.assertEqual(result, ">seq1\nAAAAAAAA\n")

    seq = self.fake_sequence(2, 'G')
    success = splitter.write_sequence_to_cluster(sequence_cluster_map, seq)
    self.assertTrue(success)
    result_file = splitter.cluster_output_files['cluster_A']
    result_file.seek(0)
    result = result_file.read()
    self.assertEqual(result, ">seq1\nAAAAAAAA\n>seq2\nGGGGGGGG\n")

    self.assertEqual(splitter.cluster_output_files['cluster_B'].read(), '')

    seq = self.fake_sequence(3, 'T')
    success = splitter.write_sequence_to_cluster(sequence_cluster_map, seq)
    self.assertTrue(success)
    result_file = splitter.cluster_output_files['cluster_B']
    result_file.seek(0)
    result = result_file.read()
    self.assertEqual(result, ">seq3\nTTTTTTTT\n")

    seq = self.fake_sequence(4, 'C')
    success = splitter.write_sequence_to_cluster(sequence_cluster_map, seq)
    # Doesn't raise an exception
    self.assertFalse(success)

  @patch('biopericles.ClusterSplitter.Bio', create=True)
  @patch('biopericles.ClusterSplitter.open', create=True)
  def test_write_all_sequences(self, open_mock, biopython_mock):

    sequences = [self.fake_sequence(*args) for args in zip(range(1,5), 'ACGT')]
    biopython_mock.SeqIO.parse.return_value = (seq for seq in sequences)

    fake_multifasta_file = StringIO()
    open_mock.return_value = fake_multifasta_file

    splitter = self.uninitialised_splitter()
    splitter.multifasta_path = '/home/multifasta.aln'
    splitter.sequence_cluster_map = {}
    splitter.write_sequence_to_cluster = MagicMock(return_value=True)

    splitter.write_all_sequences()

    open_mock.assert_called_once_with('/home/multifasta.aln', 'r')
    biopython_mock.SeqIO.parse.assert_called_once_with(fake_multifasta_file, 'fasta')
    for seq in sequences:
      splitter.write_sequence_to_cluster.assert_any_call({}, seq)

    expected_success = {'seq1': 1, 'seq2': 1, 'seq3': 1, 'seq4': 1}
    self.assertEqual(splitter.sequence_write_success, expected_success)

  @patch('biopericles.ClusterSplitter.Bio', create=True)
  @patch('biopericles.ClusterSplitter.open', create=True)
  def test_write_all_sequences_duplicates_and_missing(self, open_mock, biopython_mock):

    sequences = [self.fake_sequence(*args) for args in zip(range(1,5), 'ACGT')]
    sequences.append(sequences[0])
    biopython_mock.SeqIO.parse.return_value = (seq for seq in sequences)

    fake_multifasta_file = StringIO()
    open_mock.return_value = fake_multifasta_file

    splitter = self.uninitialised_splitter()
    splitter.multifasta_path = '/home/multifasta.aln'
    sequence_cluster_map = {
                                         'seq1': 'cluster_A',
                                         'seq2': 'cluster_A',
                                         'seq3': 'cluster_A',
                                         'seq4': 'cluster_A',
                                         'seq_another': 'cluster_B'
                                       }
    splitter.sequence_cluster_map = sequence_cluster_map
    splitter.write_sequence_to_cluster = MagicMock(return_value=True)

    splitter.write_all_sequences()

    open_mock.assert_called_once_with('/home/multifasta.aln', 'r')
    biopython_mock.SeqIO.parse.assert_called_once_with(fake_multifasta_file, 'fasta')
    for seq in sequences:
      splitter.write_sequence_to_cluster.assert_any_call(sequence_cluster_map, seq)

    expected_success = {'seq1': 2, 'seq2': 1, 'seq3': 1, 'seq4': 1, 'seq_another': 0}
    self.assertEqual(splitter.sequence_write_success, expected_success)
