import os

class BadMetadataException(ValueError):
  pass

class SampleMetadataReader(object):

  def __init__(self, filename):
    if not self.config_file_exists(filename):
      raise ValueError("Could not find filename '%s'" % filename)

  def config_file_exists(self, filename):
    return os.path.isfile(filename)

  def parse_line(self, line_data, sample_name_idx, cluster_name_idx):
    return (line_data[sample_name_idx], line_data[cluster_name_idx])

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
