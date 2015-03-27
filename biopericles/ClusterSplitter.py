import os

class NotDirectoryException(ValueError):
  pass

class ClusterSplitter(object):
  def __init__(self, multifasta, sequence_to_cluster_map, output_directory=None):
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
