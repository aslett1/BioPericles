import os
import random
import logging

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
    set_per_base = self._get_set_per_base(sequence_strings)
    bases = []
    for i,nucleotide_set in enumerate(set_per_base):
      self.logger.info("Calculating the base in the %sth position" % i)
      bases.append(self._compare_nucleotides(nucleotide_set))
    if len(bases) != self._get_length_of_longest(sequence_strings):
      raise ValueError("One sequence longer than the others, please make sure they are aligned")
    return "".join(bases)

  def _compare_nucleotides(self, neucleotide_set):
    try:
      (nucleotide,) = neucleotide_set
      if nucleotide == '-':
        return 'N'
      else:
        return nucleotide
    except ValueError:
      # More than one base in this position
      return 'N'

  def _get_all_sequences(self, sequence_generator):
    return [sequence.seq for sequence in sequence_generator]

  def _get_set_per_base(self, sequence_strings):
    return [set(bases) for bases in zip(*sequence_strings)]

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
