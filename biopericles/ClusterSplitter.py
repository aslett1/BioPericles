import Bio
import os

from collections import Counter

class NotDirectoryException(ValueError):
  pass

class ClusterSplitter(object):
  def __init__(self, multifasta, sequence_to_cluster_map, output_directory=None):
    self.cluster_output_files = {}

    if output_directory == None:
      output_directory = os.getcwd()
    self.output_directory = self.absolute_directory_path(output_directory)

  def absolute_directory_path(self, path):
    if not os.path.isdir(path):
      raise NotDirectoryException("output_directory '%s' is not a directory" % path)
    directory_name = os.path.normpath(path) # removes rogue '/' from path
    return os.path.abspath(directory_name)

  def get_clusters(self, sequence_to_cluster_map):
    clusters = sequence_to_cluster_map.values()
    return sorted(list(set(clusters)))

  def get_sequences(self, sequence_to_cluster_map):
    sequences = sequence_to_cluster_map.keys()
    return sorted(list(set(sequences)))

  def create_cluster_output_files(self, clusters):
    def create_file(cluster):
      filename = "cluster_{cluster}_multifasta.aln".format(cluster=cluster)
      path = os.path.join(self.output_directory, filename)
      return open(path, 'w')
    try:
      self.cluster_output_files = { cluster: create_file(cluster) for cluster in clusters }
    except IOError as e:
      self.output_files = {}
      raise IOError("Could not open output file to write cluster.  Exception was:\n%s" % e)

  def write_sequence_to_cluster(self, sequence_to_cluster_map, seq):
    sequence_name = seq.id
    output_cluster = sequence_to_cluster_map.get(sequence_name)
    output_file = self.cluster_output_files.get(output_cluster)
    if output_file is not None:
      output_file.write(seq.format('fasta'))
      return True
    return False

  def write_all_sequences(self):
    names_of_sequences = self.get_sequences(self.sequence_to_cluster_map)
    self.sequence_write_success = Counter({sequence: 0 for sequence in names_of_sequences})
    multifasta_file = open(self.multifasta_path, 'r')
    sequences = Bio.SeqIO.parse(multifasta_file, 'fasta')
    for seq in sequences:
      success = self.write_sequence_to_cluster(self.sequence_to_cluster_map, seq)
      if success:
        self.sequence_write_success[seq.id] += 1
    return dict(self.sequence_write_success)
