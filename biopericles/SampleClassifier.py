class SampleClassifier(object):
  def classify(self, features, feature_labels):
    """Takes a sample and classifies it"""
    pass

  def export(self, output_file):
    """Writes the classifier to a file"""
    pass

  def load(self, input_file):
    """Loads a classifier from a file"""
    pass

class BuildSampleClassifier(object):
  def __init__(self):
    self.sample_to_cluster_map = None
    self.all_features_file = None
    self.relevant_feature_labels = []
    self.training_data = None
    self.test_data = None

  def train(self):
    """Trains a classifier using the training data"""
    pass

  def get_classifier(self):
    """Returns a SampleClassifier which can be used elsewhere"""
    pass

  def test(self):
    """Returns statistics about the current classifier's performance on test
    data"""
    pass

  def get_feature_labels_from_file(self, feature_file):
    """Parses a file of features and returns the feature labels

    Some files have more features others don't.  This is a useful
    pre-computation step which extracts the feature labels used in one file
    which can be used to filter the features in another"""
    pass

  def load_data(self, features_file, sample_to_cluster_map, relevant_features):
    """Loads features from a CSV and labels the samples by cluster

    Also removes features which are not in a specified whitelist"""
    # load the data from a file of features
    # keep a copy of the feature labels
    # remove the irrelevant features
    # split the data into training and test sets
    pass

  def _extract_features(self, features_file):
    """Parses a feature file and returns the features and feature labels"""
    pass

  def _label_data(self, features, sample_to_cluster_map):
    """Takes samples and a mapping and labels the samples"""
    pass

  def _split_data(self, data, test_split=0.3):
    """Splits data into training and test sets
    
    Uses test_split to determine what proportion of the available data should be
    held back for testing (defaults to 30%) and how much should be used for
    training"""
    pass
