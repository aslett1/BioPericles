import Bio.SeqIO

class ClusterFastaFile(object):

  def __init__(self):
    self.sequences = None

  def get_consensus(self):
    return self.calculate_consensus(self.sequences)

  def load_file(self, cluster_file):
    self.sequences = Bio.SeqIO.parse(cluster_file, 'fasta')

  def calculate_consensus(self, sequence_generator):
    consensus = sequence_generator.next()
    for sequence in sequence_generator:
      consensus = self.get_pair_consensus(consensus, sequence)
    return consensus

  def get_pair_consensus(self, sequence_1, sequence_2):
    if len(sequence_1) != len(sequence_2):
      message = """\
Cannot evaluate consensus between sequences of differing lengths.
Sequence 1 was {} characters long, Sequence 2 was {} characters long.
"""
      raise ValueError(message.format(len(sequence_1), len(sequence_2)))
    consensus = map(self.compare_nucleotides, *(sequence_1,sequence_2))
    return ''.join(consensus)

  def compare_nucleotides(self, n1, n2):
    if n1 == '-' or n2 == '-':
      return 'N'
    elif n1 == n2:
      return n1
    else:
      return 'N'
