from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import train_test_split
import numpy as np
import csv

def xor(a,b):
  return (a | b) & ((a & b)==False)

feature_labels = np.array(['round', 'green', 'red', 'thick_skin', 'pips',
                           'big'])
green_apple = [1,1,0,0,1,0]
red_apple = [1,0,1,0,1,0]
orange = [1,0,0,1,1,1]
seedless_orange = [1,0,0,1,0,1]
pear = [0,1,0,0,1,0]
banana = [0,0,0,1,0,1]

perfect_features = np.array([green_apple]*100 + [red_apple]*100 + [orange]*100 + \
                    [seedless_orange]*100 + [pear]*200 + [banana]*200)
labels = np.array(['apple']*200 + ['orange']*200 + ['pear']*200 + \
                  ['banana']*200)

fuzz = np.random.rand(*perfect_features.shape) < 0.01
features = xor(perfect_features, fuzz)

sample_names = np.array([["sample_%s" % i for i in range(features.shape[0])]]).T
header_row = np.concatenate(([['Features']], [feature_labels]), axis=1)

data = np.concatenate((header_row, np.concatenate((sample_names, features), axis=1)), axis=0)
metadata = np.concatenate((sample_names, np.array([labels]).T), axis=1)
metadata = np.concatenate(([['Sample', 'Label']], metadata), axis=0)

np.savetxt('features.csv', data, fmt='%s', delimiter=',')
np.savetxt('metadata.csv', metadata, fmt='%s', delimiter=',')

train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size=0.3)

rf = RandomForestClassifier()
rf.fit(train_features, train_labels)
print rf.score(test_features, test_labels)
