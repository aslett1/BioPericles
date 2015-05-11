import unittest

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
