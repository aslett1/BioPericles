import os

class ClusterSplitter(object):

  def create_default_output_directory(self, multifasta):
    return os.path.dirname(os.path.abspath(multifasta))
