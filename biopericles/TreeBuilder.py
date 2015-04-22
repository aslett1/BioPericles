import Bio.SeqIO
import Bio.Phylo
import os
import re
import shutil
import tempfile

from biopericles.Common import LoadFastaMixin, \
                               ExternalApplicationException, \
                               RunExternalApplicationMixin

class RaxmlException(ExternalApplicationException):
  pass

class FastmlException(ExternalApplicationException):
  pass

class TreeBuilder(LoadFastaMixin, RunExternalApplicationMixin):
  def __init__(self):
    self.sequences = None # a dictionary of {sequence_name: SeqIO object}
    self.tree = None # a BioPython Phylo tree
    self.sequences_output_file = None
    self.tree_output_file = None

  def build_tree(self):
    phylip = self._create_temporary_phylip(self.sequences)
    phylip.close()
    output_directory = tempfile.mkdtemp()
    raxml_stdout, raxml_stderr = self._run_raxml('raxmlHPC', {}, phylip.name,
                                                 output_directory)
    raxml_tree_filename = self._get_raxml_tree_file(raxml_stdout)
    (tree,) = Bio.Phylo.parse(raxml_tree_filename, 'newick')
    tree.root_at_midpoint()
    self.tree = tree
    os.remove(phylip.name)
    shutil.rmtree(output_directory)
    return tree

  def add_hereditary_nodes(self):
    """Adds hereditary nodes to self.tree and self.sequences"""
    output_directory = tempfile.mkdtemp()

    temporary_sequence_file = tempfile.NamedTemporaryFile(dir=output_directory,
                                                         delete=False)
    temporary_tree_file = tempfile.NamedTemporaryFile(dir=output_directory,
                                                      delete=False)

    self._write_sequences(self.sequences.values(), temporary_sequence_file)
    self._write_tree(self.tree, temporary_tree_file)

    temporary_sequence_file.close()
    temporary_tree_file.close()

    self._run_fastml('fastml', {},
                     temporary_tree_file.name,
                     temporary_sequence_file.name,
                     output_directory)

    output_sequences_filename = os.path.join(output_directory,
                                             'all_nodes.mfa')
    output_tree_filename = os.path.join(output_directory,
                                        'all_nodes.newick')

    with open(output_sequences_filename, 'r') as output_sequences_file:
      self.load_fasta_sequences(output_sequences_file)
    self.tree = Bio.Phylo.read(output_tree_filename, 'newick')

    os.remove(temporary_sequence_file.name)
    os.remove(temporary_tree_file.name)
    shutil.rmtree(output_directory)

  def write_output(self):
    self._write_tree(self.tree, self.tree_output_file)
    self._write_sequences(self.sequences.values(), self.sequences_output_file)

  def _write_tree(self, tree, output_file):
    Bio.Phylo.write(tree, output_file, 'newick')

  def _write_sequences(self, sequences, output_file):
    Bio.SeqIO.write(sequences, output_file, 'fasta')

  def _create_temporary_phylip(self, sequences):
    """Outputs sequence dictionary as a phylip and returns a filehandle"""
    output_file = tempfile.NamedTemporaryFile('w', delete=False)
    Bio.SeqIO.write(self.sequences.values(), output_file, 'phylip')
    return output_file

  def _run_raxml(self, raxml_executable, raxml_arguments, phylip_filename, output_directory):
    """Returns the stdout from raxml, raises an exception if it exits badly"""
    default_arguments = {
      "-m": "GTRGAMMA",
      "-p": "12345",
      "-s": phylip_filename,
      "-n": "raxml_output",
      "-w": os.path.abspath(output_directory)
    }
    stdout, stderr, returncode = self._run_application(raxml_executable,
                                                       default_arguments,
                                                       raxml_arguments)
    if returncode != 0:
      raise RaxmlException("Problem running raxml on %s; some output in %s" %
                         (phylip_filename, output_directory),
                           returncode,
                           stdout,
                           stderr)
    return (stdout, stderr)

  def _get_raxml_tree_file(self, raxml_stdout):
    """Gets the filename for the raxml generated tree from the raxml output"""
    try:
      (match,) = re.finditer("Best-scoring ML tree written to:\s+(\S+)",
                             raxml_stdout)
      return match.group(1) # the path to the output file
    except ValueError:
      # There wasn't just one match so we can't trust the output
      return None

  def _run_fastml(self, fastml_executable, fastml_arguments, tree_filename, sequence_fasta_filename, output_directory):
    default_arguments = {
      "-s": sequence_fasta_filename,
      "-t": tree_filename,
      "-x": os.path.join(output_directory, 'all_nodes.newick'),
      "-j": os.path.join(output_directory, 'all_nodes.mfa'),
      "-mg": '',
      "-qf": ''
    }

    stdout, stderr, returncode = self._run_application(fastml_executable,
                                                       default_arguments,
                                                       fastml_arguments,
                                                       cwd=output_directory)
    if returncode != 0:
      raise FastmlException("Problem running fastml using %s and %s; some output in %s" %
                           (sequence_fasta_filename, tree_filename, output_directory),
                           returncode,
                           stdout,
                           stderr)
    return (stdout, stderr)
