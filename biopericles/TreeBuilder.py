import Bio.SeqIO
import os
import tempfile

def try_and_get_filename(filehandle):
  try:
    return filehandle.name
  except AttributeError:
    return "ANONYMOUS_FILE"

class TreeBuilder(object):
  def __init__(self):
    self.sequences = None # a dictionary of {sequence_name: SeqIO object}

  def load_fasta_sequences(self, fasta_file):
    """Load sequences from a fasta_file into a dictionary of {sequence_name: SeqIO object}

    Takes a filehandle to a aligned multifasta
    """
    list_of_sequences = Bio.SeqIO.parse(fasta_file, 'fasta')
    self.sequences = {}
    for seq in list_of_sequences:
      if seq.name in self.sequences:
        message = "Tried to load sequences from {filename} but got multiple sequences for {sequence}"
        raise ValueError(message.format(filename=try_and_get_filename(fasta_file),
                                        sequence=seq.name))
      self.sequences[seq.name] = seq

  def build_tree(self):
    pass # returns a Bio Python tree object

  def add_hereditary_nodes(self, tree, sequences):
    pass # takes a tree, a dictionary of BioPython sequences; returns a tree object

  def write_output(self, output_directory, filename_prefix):
    pass # writes the sequences and internal nodes and full tree

  def _create_temporary_phylip(self, fasta_file):
    pass # reformats the fasta file as a phylip and returns a filehandle

  def _run_raxml(self, phylip_filename):
    pass # returns the stdout from raxml, raises an exception if it exits badly

  def _get_tree_filename(self, raxml_stdout):
    pass # gets the filename for the tree from the raxml output

  def _run_fastml(self, output_directory, tree_filename, sequence_fasta_filename):
    pass # runs fastml using the tree from raxml and the sequence

  def _find_fastml_tree_file(self, output_directory):
    pass # return a filehandle for the fastml output tree

  def _find_fastml_sequence_file(self, output_directory):
    tpass # return a filehandle for the fastml output sequences
