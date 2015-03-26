import os

class SampleMetadataReader(object):

  def __init__(self, filename):
    if not self.config_file_exists(filename):
      raise ValueError("Could not find filename '%s'" % filename)

  def config_file_exists(self, filename):
    return os.path.isfile(filename)
