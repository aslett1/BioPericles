import os
import numpy as np
import unittest

from mock import MagicMock, patch
from numpy.testing import assert_array_equal
from sklearn.ensemble import RandomForestClassifier
from StringIO import StringIO

from biopericles.SampleClassifier import SampleClassifier, \
                                         BuildSampleClassifier, \
                                         SortFeaturesMixin
from biopericles.Common import context_aware_tempfile

def test_data():
  this_file = os.path.abspath(__file__)
  this_dir = os.path.dirname(this_file)
  return os.path.join(this_dir, 'data')

class TestSortFeaturesMixin(unittest.TestCase):
  def test_sort_features(self):
    sorter = SortFeaturesMixin()

    features = np.array([[0,0,0,1], [0,0,1,0], [0,0,1,1], [0,1,0,0]])
    old_labels = np.array(['feature_1', 'feature_2', 'feature_3', 'feature_4'])
    new_labels = np.array(['feature_1', 'feature_3', 'feature_2',
                           'unknown_feature'])

    new_features, labels, missing_labels = sorter.sort_features(features,
                                                                old_labels,
                                                                new_labels)

    expected_new_features = np.array([[0,0,0], [0,1,0], [0,1,0], [0,0,1]])
    expected_labels = np.array(['feature_1', 'feature_3', 'feature_2'])
    expected_missing_labels = np.array(['unknown_feature'])

    assert_array_equal(new_features, expected_new_features)
    assert_array_equal(labels, expected_labels)
    assert_array_equal(missing_labels, expected_missing_labels)

class TestSampleClassifier(unittest.TestCase):
  def check_latest(self, mock, expected_array):
    ((arg,), kwargs) = mock.call_args
    assert_array_equal(arg, expected_array)

  def test_classify(self):
    rf = RandomForestClassifier()
    rf.fit(np.array([[0,1]]), np.array(['cluster_1']))

    classifier_mock = MagicMock()
    classifier_mock.predict.side_effect = rf.predict
    classifier = SampleClassifier(classifier_mock, np.array(['feature_1',
                                                             'feature_2']))
    classifier.logger = MagicMock()

    classifier.classify(np.array([0,1]))
    self.check_latest(classifier.classifier.predict, np.array([[0,1]]))

    classifier.classify(np.array([[0,1],[1,0]]))
    self.check_latest(classifier.classifier.predict, np.array([[0,1],[1,0]]))

    labels = np.array(['feature_2', 'feature_1'])
    classifier.classify(np.array([0,1]), feature_labels=labels)
    self.check_latest(classifier.classifier.predict, np.array([[1,0]]))

    labels = np.array(['feature_1', 'feature_2'])
    classifier.classify(np.array([0,1]), feature_labels=labels)
    self.check_latest(classifier.classifier.predict, np.array([[0,1]]))

    labels = np.array(['feature_1', 'unknown_label'])
    self.assertRaises(ValueError, classifier.classify, np.array([0,1]), feature_labels=labels)
    classifier.logger.error.assert_any_call("Could not find feature 'feature_2', ignoring it")

    self.assertRaises(ValueError, classifier.classify, np.array([0]))
    classifier.logger.error.assert_any_call("Classifier trained with 2 features, got 1 to make prediction with")
    self.assertRaises(ValueError, classifier.classify, np.array([0,0,1]))
    classifier.logger.error.assert_any_call("Classifier trained with 2 features, got 3 to make prediction with")

    labels = np.array(['feature_1'])
    self.assertRaises(ValueError, classifier.classify, np.array([0,1]), feature_labels=labels)
    classifier.logger.error.assert_any_call("Classifier trained with 2 features, got 1 to make prediction with")

  @patch('biopericles.SampleClassifier.logging')
  def test_load(self, logging_mock):
    # I've patched logging because I don't want the noise
    path_to_classifier = os.path.join(test_data(), 'classifier.pkl')
    actual_md5sum = '9c829ce1a9806a6975442e761f4e616d'
    actual_sha256sum = 'b7e9ad6835fd3f4970229c4c5ab97a14bd2c118da3f2b73ce9df10072b7aa597'
    with open(path_to_classifier, 'rb') as classifier_file:
      self.assertRaises(ValueError, SampleClassifier.load, classifier_file)
      self.assertRaises(ValueError, SampleClassifier.load, classifier_file,
                        md5sum='wrong')
      self.assertRaises(ValueError, SampleClassifier.load, classifier_file,
                        sha256sum='wrong')
      self.assertRaises(ValueError, SampleClassifier.load, classifier_file,
                        md5sum='wrong', sha256sum=actual_sha256sum)
      classifier = SampleClassifier.load(classifier_file, md5sum=actual_md5sum)
      self.assertIsInstance(classifier.classifier, RandomForestClassifier)
      classifier = SampleClassifier.load(classifier_file, sha256sum=actual_sha256sum)
      classifier = SampleClassifier.load(classifier_file, force=True)
      classifier = SampleClassifier.load(classifier_file, md5sum='wrong', force=True)
      assert_array_equal(classifier.feature_labels, np.array(['feature_1',
                                                              'feature_2']))

  def test_export(self):
    rf = RandomForestClassifier()
    rf.fit(np.array([[0,1]]), np.array(['cluster_1']))

    classifier = SampleClassifier(rf, np.array(['feature_1','feature_2']))
    classifier.logger = MagicMock()

    with context_aware_tempfile('w', delete=False) as output_file:
      classifier.export(output_file)
      output_file.close()
      input_file = open(output_file.name, 'r')
      new_classifier = SampleClassifier.load(input_file, force=True)

################################################################################
##   Uncomment me if you want to create a new classifier test file.  Remember to
##   update the md5sum and sha256sum values in the test_load function above
##   afterwards even if you haven't changed how the classifier is created.
################################################################################
#
#    path_to_classifier = os.path.join(test_data(), 'classifier.pkl')
#    with open(path_to_classifier, 'wb') as output_file:
#      classifier.export(output_file)
#
################################################################################

    self.assertIsInstance(new_classifier.classifier, RandomForestClassifier)
    assert_array_equal(new_classifier.feature_labels, np.array(['feature_1',
                                                                'feature_2']))

class TestBuildSampleClassifier(unittest.TestCase):
  @patch('biopericles.SampleClassifier.np.random')
  def test_load_data(self, random_mock):
    feature_builder = BuildSampleClassifier()
    feature_builder.logger = MagicMock()

    file_contents = """\
Features,feature_1,feature_2,feature_3
sample_1,0,1,1
sample_2,1,0,1
sample_3,1,0,0
"""
    feature_file = StringIO(file_contents)
    sample_to_cluster_map = {'sample_1': 'cluster_A', 'sample_2': 'cluster_B',
                             'unused_sample': 'unused_cluster'}
    relevant_features = np.array(['feature_1', 'feature_2', 'missing_feature'])

    random_mock.uniform.return_value = np.array([0.1,0.9])

    feature_builder.load_data(feature_file, sample_to_cluster_map,
                              relevant_features)

    expected_training_labels = np.array(['cluster_A'])
    expected_training_features = np.array([[0, 1]])

    expected_testing_labels = np.array(['cluster_B'])
    expected_testing_features = np.array([[1, 0]])

    assert_array_equal(feature_builder.training_labels,
                       expected_training_labels)
    assert_array_equal(feature_builder.training_features,
                       expected_training_features)
    assert_array_equal(feature_builder.testing_labels,
                       expected_testing_labels)
    assert_array_equal(feature_builder.testing_features,
                       expected_testing_features)

    feature_builder.logger.warn.assert_any_call("Could not assign sample 'sample_3' to a cluster, skipping")
    feature_builder.logger.warn.assert_any_call("Could not find feature 'missing_feature', ignoring it")

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

  def test_samples_with_unknown_clusters(self):
    feature_builder = BuildSampleClassifier()

    sample_names = np.array(['sample_1', 'unknown_sample', 'sample_3'])
    cluster_labels = np.array(['cluster_A', None, 'cluster_A'])

    expected_sample_names = np.array(['unknown_sample'])
    actual_sample_names = feature_builder._samples_with_unknown_clusters(sample_names,
                                                                         cluster_labels)

    assert_array_equal(actual_sample_names, expected_sample_names)

  def test_only_labeled_data(self):
    feature_builder = BuildSampleClassifier()

    features = np.array([[0,0,1], [0,1,0], [1,0,0]])
    cluster_labels = np.array(['cluster_A', None, 'cluster_A'])

    expected_features = np.array([[0,0,1], [1,0,0]])
    expected_cluster_labels = np.array(['cluster_A', 'cluster_A'])

    actual_cluster_labels, actual_features = feature_builder._only_labeled_data(cluster_labels, features)

    assert_array_equal(actual_cluster_labels, expected_cluster_labels)
    assert_array_equal(actual_features, expected_features)

  def test_warn_about_missing_features(self):
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
