import csv
import logging
import numpy as np

from sklearn.ensemble import RandomForestClassifier

from biopericles.Common import try_and_get_filename

class SampleClassifier(object):
  def __init__(self, classifier, feature_labels):
    self.classifier = classifier
    self.feature_labels = feature_labels

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
    self.relevant_feature_labels = []
    self.feature_labels = None
    self.training_labels = None
    self.training_features = None
    self.testing_labels = None
    self.testing_features = None
    self.logger = logging.getLogger(__name__)
    self.classifier = None

  def train(self, *args, **kwargs):
    """Trains a classifier using the training data"""
    rf = RandomForrestClassifier(*args, **kwargs)
    rf.fit(self.training_features, self.training_labels)
    self.classifier = SampleClassifier(rf, self.feature_labels)

  def get_classifier(self):
    """Returns a SampleClassifier which can be used elsewhere"""
    return self.classifier

  def test(self):
    """Returns statistics about the current classifier's performance on test
    data"""
    assert self.feature_labels == self.classifier.feature_labels
    return self.classifier.classifier.score(self.testing_features, self.testing_labels)

  def get_feature_labels_from_file(self, feature_file):
    """Parses a file of features and returns the feature labels

    Some files have more features others don't.  This is a useful
    pre-computation step which extracts the feature labels used in one file
    which can be used to filter the features in another"""
    feature_file.seek(0)
    features_line = csv.reader(feature_file).next()
    if features_line[0] != 'Features':
      filename = try_and_get_filename(feature_file)
      raise ValueError("Issue parsing '%s', expected first element to be 'Features'" % filename)
    feature_labels = features_line[1:]
    return np.array(feature_labels)

  def load_data(self, features_file, sample_to_cluster_map, relevant_features):
    """Loads features from a CSV and labels the samples by cluster

    Also removes features which are not in a specified whitelist"""
    # load the data from a file of features
    sample_names, features = self._extract_features(features_file)
    # keep a copy of the feature labels
    cluster_labels = self._get_cluster_labels_for_samples(sample_names,
                                                          sample_to_cluster_map)
    # get feature labels
    feature_labels = self.get_feature_labels_from_file(features_file)
    # remove the irrelevant features
    features, feature_labels, missing_feature_labels = self._only_relevant_features(features,
                                                                                    feature_labels,
                                                                                    relevant_features)
    # warn about missing features
    self._warn_about_missing_features(missing_feature_labels)
    # remove unlabeled data (and warn the user)
    cluster_labels_for_features = self._get_cluster_labels_for_samples(sample_names,
                                                                sample_to_cluster_map)
    unlabeled_samples = self._samples_with_unknown_clusters(sample_names,
                                                            cluster_labels_for_features)
    self._warn_about_unlabeled_samples(unlabeled_samples)
    clusters_with_labels, features_with_labels = self._only_labeled_data(cluster_labels_for_features,
                                                                         features)
    # split the data into training and test sets
    split = self._split_data(clusters_with_labels, features_with_labels, test_split=0.3)
    self.training_labels, self.training_features = split[:2]
    self.testing_labels, self.testing_features = split[2:]
    self.feature_labels = feature_labels

  def _check_features_binary(self, features):
    check_zero_or_one = np.vectorize(lambda v: (v==1) or (v==0))
    return check_zero_or_one(features).all()

  def _extract_features(self, features_file):
    """Parses a feature file and returns the sample_names and features"""
    features_file.seek(0)
    filename = try_and_get_filename(features_file)
    feature_file_lines = [line for line in csv.reader(features_file)]
    data = np.array(feature_file_lines)
    try:
      number_of_rows, number_of_columns = data.shape
      # shape returns a 2 element tuple if there are a consistent number of
      # columns.  Otherwise it returns a 1 element tuple which would trigger a
      # ValueError.
    except ValueError as e:
      raise ValueError("Issue parsing '%s', some samples have more features than others" % filename)
    feature_labels = data[0,1:]
    sample_names = data[1:,0]
    str_to_int = np.vectorize(int)
    try:
      features = str_to_int(data[1:,1:])
    except ValueError:
      raise ValueError("Issue parsing '%s', expected features to be 0 or 1" % filename)
    if not self._check_features_binary(features):
      raise ValueError("Issue parsing '%s', expected features to be 0 or 1" % filename)
    return sample_names, features

  def _get_cluster_labels_for_samples(self, sample_names, sample_to_cluster_map):
    """Takes samples and a mapping and labels the samples"""
    sample_names = sample_names.copy()
    sample_to_cluster = np.vectorize(sample_to_cluster_map.get, otypes=[object])
    return sample_to_cluster(sample_names)

  def _only_relevant_features(self, features, feature_labels,
                              relevant_feature_labels):
    """Removes features which are not 'relevant'

    Returns features, feature_labels which were found in relevant_features and
    relevant_features which could not be found"""
    feature_labels_list = list(feature_labels)
    indices_of_common_features = [feature_labels_list.index(label) for label in
                                  relevant_feature_labels if label in
                                  feature_labels_list]
    common_features = features[:,indices_of_common_features]
    common_feature_labels = feature_labels[indices_of_common_features]
    is_missing_feature = np.in1d(relevant_feature_labels, feature_labels) == False
    missing_feature_labels = relevant_feature_labels[is_missing_feature]
    return common_features, common_feature_labels, missing_feature_labels

  def _warn_about_missing_features(self, features):
    for feature in features:
      self.logger.warn("Could not find feature '%s', ignoring it" % feature)

  def _only_labeled_data(self, cluster_labels, features):
    is_not_none = np.vectorize(lambda v: v is not None)
    known_labels = is_not_none(cluster_labels)
    clusters_with_labels = cluster_labels[known_labels]
    features_with_labels = features[known_labels]
    return clusters_with_labels, features_with_labels

  def _warn_about_unlabeled_samples(self, samples):
    for sample in samples:
      self.logger.warn("Could not assign sample '%s' to a cluster, skipping" %
                       sample)

  def _samples_with_unknown_clusters(self, sample_names, cluster_labels):
    is_none = np.vectorize(lambda v: v is None)
    unknown_labels = is_none(cluster_labels)
    return sample_names[unknown_labels]

  def _split_data(self, labels, features, test_split=0.3):
    """Splits data into training and test sets

    Uses test_split to determine what proportion of the available data should be
    held back for testing (defaults to 30%) and how much should be used for
    training"""

    if len(labels) != len(features):
      raise ValueError("Cannot split data without equal number of sample names and features")

    training_split = np.random.uniform(0,1,len(labels)) < (1.0 - test_split)
    testing_split = training_split == False

    training_labels = labels[training_split]
    training_features = features[training_split]

    testing_labels = labels[testing_split]
    testing_features = features[testing_split]

    return training_labels, training_features, testing_labels, testing_features
