import os
import csv

class BadMetadataException(ValueError):
  pass

class SampleMetadataReader(object):

  def __init__(self, metadata_filename, sample_name_idx, cluster_name_idx):
    self.cluster_sample_map = None
    self.sample_cluster_map = None
    if not self.config_file_exists(metadata_filename):
      raise ValueError("Could not find filename '%s'" % metadata_filename)
    self.sample_name_idx = sample_name_idx
    self.cluster_name_idx = cluster_name_idx
    with open(metadata_filename, 'r') as metadata_file:
      self.parse(metadata_file)

  def config_file_exists(self, filename):
    return os.path.isfile(filename)

  def parse(self, metadata_file):
    csv_reader = csv.reader(metadata_file, delimiter=',')
    header = csv_reader.next()
    sample_cluster_tuples = map(self.parse_line, csv_reader)
    self.cluster_sample_map = self.create_cluster_sample_map(sample_cluster_tuples)
    self.sample_cluster_map = self.create_sample_cluster_map(sample_cluster_tuples)

  def parse_line(self, line_data):
    return (line_data[self.sample_name_idx], line_data[self.cluster_name_idx])

  def create_cluster_sample_map(self, sample_cluster_tuples):
    cluster_sample_map = {}
    for sample, cluster in sample_cluster_tuples:
      cluster_sample_map.setdefault(cluster, []).append(sample)
    unique_cluster_sample_map = {k: sorted(list(set(v))) for k,v in cluster_sample_map.items()}
    return unique_cluster_sample_map

  def create_sample_cluster_map(self, sample_cluster_tuples):
    deduplicated_sample_cluster_tuples = list(set(sample_cluster_tuples))
    samples = list(set([sample for sample, cluster in deduplicated_sample_cluster_tuples]))
    if len(samples) != len(deduplicated_sample_cluster_tuples):
      raise BadMetadataException("Each sample should only be attributed to one cluster")
    return { sample: cluster for sample, cluster in deduplicated_sample_cluster_tuples }
