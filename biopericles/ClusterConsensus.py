import logging
import os
import random

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ClusterConsensus(object):

  def __init__(self):
    self.sequences = None # a SeqRecord generator or list
    self.output_file = None
    self.cluster_name = None
    self.logger = logging.getLogger(__name__)

  def load_fasta_file(self, cluster_file):
    self.sequences = Bio.SeqIO.parse(cluster_file, 'fasta')
    if not self.cluster_name:
      self.cluster_name = self._get_cluster_name()

  def get_consensus(self):
    return self._calculate_consensus(self.sequences)

  def create_consensus_sequence(self):
    consensus = self._calculate_consensus(self.sequences)
    seq = Seq(consensus)
    name = self._get_cluster_name()
    return SeqRecord(seq, description='', id=name)

  def write_consensus(self):
    consensus_sequence = self.create_consensus_sequence()
    Bio.SeqIO.write(consensus_sequence, self.output_file, 'fasta')


  def _calculate_consensus(self, sequence_generator):
    sequence_strings = self._get_all_sequences(sequence_generator)
    bases = []
    try:
      for i,bases_in_ith_position in self._iterate_bases(sequence_strings):
        # bases_in_the_ith_position is a list of each of the bases at a particular
        # position in each of the sequences (e.g. the first or second base in each
        # position)
        if i % 100000 == 0:
          self.logger.info("Calculating the base in the %sth position" % i)
        bases.append(self._compare_bases(bases_in_ith_position))
    except IndexError:
      # self._iterate_bases assumes that all sequences are the same length as
      # the first one.  If not the sequences cannot have been aligned properly
      # so raise an error.
      raise ValueError("One sequence longer than the others, please make sure they are aligned")
    return "".join(bases)

  def _iterate_bases(self, sequence_strings):
    """Returns a generator of (i, bases)

    For each sequence, this returns an array of the ith base.

    Given:
    sequence_strings = ['ABC', 'DEF']

    Returns a generator for:
    (0, ['A', 'D']), (1, ['B', 'E']), (2, ['C', 'F'])"""
    length_of_sequences = len(sequence_strings[0])
    return ((i, [sequence[i] for sequence in sequence_strings]) for i in xrange(length_of_sequences))

  def _compare_bases(self, bases):
    first_base = bases[0]
    if not first_base in 'AGCT':
      return 'N'
    elif len(bases) == bases.count(first_base):
      # Previously I created a set of all the bases,
      # this is 10 times quicker
      return first_base
    else:
      return 'N'

  def _get_all_sequences(self, sequence_generator):
    self.logger.info("Loading all the sequences into memory")
    return [sequence.seq for sequence in sequence_generator]

  def _get_cluster_name(self):
    if self.cluster_name:
      return self.cluster_name
    try:
      filename = self.output_file.name
      name = os.path.splitext(os.path.basename(filename))[0]
    except AttributeError:
      name = "cluster_%s" % random.randrange(10000,100000)
    return name

  def _get_length_of_longest(self, sequence_strings):
    longest = max(sequence_strings, key=len)
    return len(longest)
