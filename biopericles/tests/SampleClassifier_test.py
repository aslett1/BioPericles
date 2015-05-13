import numpy as np
import unittest

from mock import MagicMock, patch
from numpy.testing import assert_array_equal
from StringIO import StringIO

from biopericles.SampleClassifier import BuildSampleClassifier

class TestBuildSampleClassifier(unittest.TestCase):
  def test_get_feature_labels_from_file(self):
    feature_builder = BuildSampleClassifier()

    file_contents = """\
Features,feature_1,feature_2,feature_3
sample_1,0,1,1
sample_2,1,0,1
sample_3,1,0,0
"""
    feature_file = StringIO(file_contents)

    expected = ['feature_1','feature_2','feature_3']
    feature_labels = feature_builder.get_feature_labels_from_file(feature_file)
    self.assertItemsEqual(feature_labels, expected)

    # Check it seeks to the beginning of the file
    feature_labels = feature_builder.get_feature_labels_from_file(feature_file)
    self.assertItemsEqual(feature_labels, expected)

    file_contents = """\
Something,else
not,what,we,want
"""
    feature_file = StringIO(file_contents)

    self.assertRaises(ValueError, feature_builder.get_feature_labels_from_file, feature_file)

  def test_extract_features(self):
    feature_builder = BuildSampleClassifier()

    file_contents = """\
Features,feature_1,feature_2,feature_3
sample_1,0,1,1
sample_2,1,0,1
sample_3,1,0,0
"""
    feature_file = StringIO(file_contents)

    expected_samples = np.array(['sample_1','sample_2','sample_3'])
    expected_features = np.array([
      [0,1,1],
      [1,0,1],
      [1,0,0]
    ])
    samples, features = feature_builder._extract_features(feature_file)
    assert_array_equal(samples, expected_samples)
    assert_array_equal(features, expected_features)

    # Check it seeks to the beginning of the file
    samples, features = feature_builder._extract_features(feature_file)
    assert_array_equal(samples, expected_samples)
    assert_array_equal(features, expected_features)

    file_contents = """\
Features,Something,else
not,x
what,we,want
"""
    feature_file = StringIO(file_contents)
    self.assertRaises(ValueError, feature_builder._extract_features, feature_file)

    file_contents = """\
Features,feature_1,feature_2
sample_1,1,0
sample_2,1,error
"""
    feature_file = StringIO(file_contents)
    self.assertRaises(ValueError, feature_builder._extract_features, feature_file)

    file_contents = """\
Features,feature_1,feature_2
sample_1,1,2
"""
    feature_file = StringIO(file_contents)
    self.assertRaises(ValueError, feature_builder._extract_features, feature_file)

  def test_get_cluster_labels_for_samples(self):
    feature_builder = BuildSampleClassifier()

    sample_to_cluster_map = {
      'sample_1': 'cluster_A',
      'sample_2': 'cluster_A',
      'sample_3': 'cluster_A',
      'sample_4': 'cluster_B',
    }

    sample_names = np.array(['sample_1', 'sample_4', 'sample_3'])
    expected = np.array(['cluster_A', 'cluster_B', 'cluster_A'])

    clusters = feature_builder._get_cluster_labels_for_samples(sample_names,
                                                               sample_to_cluster_map)
    assert_array_equal(clusters, expected)

    sample_names = np.array(['sample_1', 'unknown_sample', 'sample_3'])
    expected = np.array(['cluster_A', None, 'cluster_A'])

    clusters = feature_builder._get_cluster_labels_for_samples(sample_names,
                                                               sample_to_cluster_map)
    assert_array_equal(clusters, expected)

  def test_only_relevant_features(self):
    feature_builder = BuildSampleClassifier()

    features = np.array([[0,0,1], [0,1,0], [1,0,0]])
    feature_labels = np.array(['feature_1', 'feature_2', 'feature_3'])

    relevant_feature_labels = np.array(['feature_1', 'feature_2'])
    expected_features = np.array([[0,0], [0,1], [1,0]])
    expected_labels = np.array(['feature_1', 'feature_2'])

    results = feature_builder._only_relevant_features(features, feature_labels,
                                                      relevant_feature_labels)
    relevant_features, labels, missing_labels = results

    assert_array_equal(relevant_features, expected_features)
    assert_array_equal(labels, expected_labels)
    self.assertEqual(list(missing_labels), [])


    features = np.array([[0,0,1], [0,1,0], [1,0,0]])
    feature_labels = np.array(['feature_1', 'feature_2', 'feature_3'])

    relevant_feature_labels = np.array(['feature_2', 'feature_1'])
    expected_features = np.array([[0,0], [1,0], [0,1]])
    expected_labels = np.array(['feature_2', 'feature_1'])
    expected_missing_labels = np.array([])

    results = feature_builder._only_relevant_features(features, feature_labels,
                                                      relevant_feature_labels)
    relevant_features, labels, missing_labels = results

    assert_array_equal(relevant_features, expected_features)
    assert_array_equal(labels, expected_labels)
    assert_array_equal(list(missing_labels), [])


    relevant_feature_labels = np.array(['feature_1', 'feature_4'])
    expected_features = np.array([[0], [0], [1]])
    expected_labels = np.array(['feature_1'])
    expected_missing_labels = np.array(['feature_4'])

    results = feature_builder._only_relevant_features(features, feature_labels,
                                                      relevant_feature_labels)
    relevant_features, labels, missing_labels = results

    assert_array_equal(relevant_features, expected_features)
    assert_array_equal(labels, expected_labels)
    assert_array_equal(missing_labels, expected_missing_labels)

  def test_only_labeled_data(self):
    feature_builder = BuildSampleClassifier()

    sample_names = np.array(['sample_1', 'unknown_sample', 'sample_3'])
    cluster_labels = np.array(['cluster_A', None, 'cluster_A'])

    expected_sample_names = np.array(['unknown_sample'])
    actual_sample_names = feature_builder._samples_with_unknown_clusters(sample_names,
                                                                         cluster_labels)

    assert_array_equal(actual_sample_names, expected_sample_names)

  def test_samples_with_unknown_clusters(self):
    feature_builder = BuildSampleClassifier()

    features = np.array([[0,0,1], [0,1,0], [1,0,0]])
    cluster_labels = np.array(['cluster_A', None, 'cluster_A'])

    expected_features = np.array([[0,0,1], [1,0,0]])
    expected_cluster_labels = np.array(['cluster_A', 'cluster_A'])

    actual_cluster_labels, actual_features = feature_builder._only_labeled_data(cluster_labels, features)

    assert_array_equal(actual_cluster_labels, expected_cluster_labels)
    assert_array_equal(actual_features, expected_features)

  def test_warn_missing_features(self):
    feature_builder = BuildSampleClassifier()
    feature_names = np.array(['unknown_feature'])

    feature_builder.logger = MagicMock()
    feature_builder._warn_about_missing_features(feature_names)

    expected_warning = "Could not find feature 'unknown_feature', ignoring it"
    feature_builder.logger.warn.assert_called_with(expected_warning)

  def test_warn_about_unlabeled_samples(self):
    feature_builder = BuildSampleClassifier()
    sample_names = np.array(['unknown_sample'])

    feature_builder.logger = MagicMock()
    feature_builder._warn_about_unlabeled_samples(sample_names)

    expected_warning = "Could not assign sample 'unknown_sample' to a cluster, skipping"
    feature_builder.logger.warn.assert_called_with(expected_warning)

  @patch('biopericles.SampleClassifier.np.random')
  def test_split_data(self, random_mock):
    feature_builder = BuildSampleClassifier()

    labels = np.array(['cluster_A', 'cluster_B', 'cluster_A', 'cluster_A', 'cluster_B'], dtype=object)
    features = np.array([[0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1]])

    random_mock.uniform.return_value = np.array([0.1,0.2,0.65,0.7,0.75])

    expected_training_labels = np.array(['cluster_A', 'cluster_B', 'cluster_A'], dtype=object)
    expected_testing_labels = np.array(['cluster_A', 'cluster_B'], dtype=object)

    expected_training_features = np.array([[0,0,1], [0,1,0], [0,1,1]])
    expected_testing_features = np.array([[1,0,0], [1,0,1]])

    results = feature_builder._split_data(labels, features, test_split=0.3)
    training_labels, training_features = results[:2]
    testing_labels, testing_features = results[2:]

    assert_array_equal(training_labels, expected_training_labels)
    assert_array_equal(training_features, expected_training_features)
    assert_array_equal(testing_labels, expected_testing_labels)
    assert_array_equal(testing_features, expected_testing_features)

    labels = np.array(['cluster_A', 'cluster_B', 'cluster_A', 'cluster_A', 'cluster_B'], dtype=object)
    features = np.array([[0,1,0], [0,1,1], [1,0,0], [1,0,1]]) # Not enough features

    self.assertRaises(ValueError, feature_builder._split_data, labels, features, test_split=0.3)
