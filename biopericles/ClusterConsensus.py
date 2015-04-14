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
    consensus = sequence_generator.next()
    for sequence in sequence_generator:
      consensus = self._get_pair_consensus(consensus, sequence)
    return consensus

  def _get_pair_consensus(self, sequence_1, sequence_2):
    if len(sequence_1) != len(sequence_2):
      message = """\
Cannot evaluate consensus between sequences of differing lengths.
Sequence 1 was {} characters long, Sequence 2 was {} characters long.
"""
      raise ValueError(message.format(len(sequence_1), len(sequence_2)))
    consensus = map(self._compare_nucleotides, *(sequence_1,sequence_2))
    return ''.join(consensus)

  def _compare_nucleotides(self, n1, n2):
    if n1 == '-' or n2 == '-':
      return 'N'
    elif n1 == n2:
      return n1
    else:
      return 'N'

  def _get_cluster_name(self):
    if self.cluster_name:
      return self.cluster_name
    try:
      filename = self.output_file.name
      name = os.path.splitext(os.path.basename(filename))[0]
    except AttributeError:
      name = "cluster_%s" % random.randrange(10000,100000)
    return name
