import os

class SampleMetadataReader(object):

  def __init__(self, filename):
    if not self.config_file_exists(filename):
      raise ValueError("Could not find filename '%s'" % filename)

  def config_file_exists(self, filename):
    return os.path.isfile(filename)

  def parse_line(self, line_data, sample_name_idx, cluster_name_idx):
    return (line_data[sample_name_idx], line_data[cluster_name_idx])
