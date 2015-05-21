#!/usr/bin/env python

import csv
import numpy as np

from biopericles.SampleClassifier import SampleClassifier

with open('/lustre/scratch108/pathogen/bt5/pericles/output/classifier.pkl', 'r') as classifier_file:
  classifier = SampleClassifier.load(classifier_file, force=True)

with open('/lustre/scratch108/pathogen/bt5/pericles/clusters.csv', 'r') as metadata_file:
  metadata = np.array([line for line in csv.reader(metadata_file)])

sample_cluster_map = {}
for row in metadata:
  sample_name = row[1]
  cluster = row[2]
  sample_cluster_map[sample_name] = cluster

with open('/lustre/scratch108/pathogen/bt5/pericles/output/all-snp-features.csv', 'r') as snp_file:
  data = np.array([line for line in csv.reader(snp_file)])

feature_labels = data[0,1:]

right, wrong = (0,0)

predicted_clusters = classifier.classify(data[1:,1:], feature_labels)

for (sample_name, predicted_cluster) in zip(data[1:,0], predicted_clusters):
  try:
    expected_cluster = sample_cluster_map[sample_name]
  except:
    print "Couldn't look up correct cluster for %s, skipping" % sample_name
    continue
  if predicted_cluster != expected_cluster:
    print "Predicted %s in %s but should be in %s" % (sample_name, predicted_cluster, expected_cluster)
    wrong += 1
  else:
    right += 1

print "Got %s wrong out of %s" % (wrong, (right+wrong))
