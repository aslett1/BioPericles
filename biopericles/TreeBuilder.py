import Bio.SeqIO
import Bio.Phylo
import logging
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
    self.logger = logging.getLogger(__name__)

  def build_tree(self):
    phylip = self._create_temporary_phylip(self.sequences)
    self.logger.info("Created temporary phylip %s for raxml" % phylip.name)
    phylip.close()
    output_directory = tempfile.mkdtemp()
    self.logger.info("Created output directory %s for raxml" % output_directory)
    raxml_stdout, raxml_stderr = self._run_raxml('raxmlHPC', {}, phylip.name,
                                                 output_directory)
    raxml_tree_filename = self._get_raxml_tree_file(raxml_stdout)
    self.logger.info("Wrote raxml tree file to %s" % raxml_tree_filename)
    if raxml_tree_filename == None:
      error_message = """\
Could not parse where raxml stored the output tree.
Raxml stdout:
%s
Raxml stderr:
%s
""" % (raxml_stdout, raxml_stderr)
      self.logger.error(error_message)
      raise RaxmlException("Could not parse where raxml stored output tree",
                           0, raxml_stdout, raxml_stderr)
    self.logger.info("raxml output tree to %s" % raxml_tree_filename)
    (tree,) = Bio.Phylo.parse(raxml_tree_filename, 'newick')
    tree.root_at_midpoint()
    self.tree = tree
    os.remove(phylip.name)
    shutil.rmtree(output_directory)
    return tree

  def add_hereditary_nodes(self):
    """Adds hereditary nodes to self.tree and self.sequences"""
    output_directory = tempfile.mkdtemp()
    self.logger.info("Created output directory %s for fastml" % output_directory)

    temporary_sequence_file = tempfile.NamedTemporaryFile(dir=output_directory,
                                                         delete=False)
    temporary_tree_file = tempfile.NamedTemporaryFile(dir=output_directory,
                                                      delete=False)

    self._write_sequences(self.sequences.values(), temporary_sequence_file)
    self.logger.info("Temporarily stored sequences in %s" %
                      temporary_sequence_file.name)
    self._write_tree(self.tree, temporary_tree_file)
    self.logger.info("Temporarily stored tree in %s" %
                      temporary_tree_file.name)

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
      self.logger.info("Reloaded ancestral sequences from %s" %
                       output_sequences_filename)
    self.tree = Bio.Phylo.read(output_tree_filename, 'newick')
    self.logger.info("Loaded ancestral tree from %s" %
                     output_tree_filename)

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
      error_message = """\
Problem running raxml on %s.
stdout:
%s
stderr:
%s
""" % (phylip_filename, stdout,stderr)
      self.logger.error(error_message)
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
    except ValueError:
      # Didn't get just one match, maybe this version of raxml is different
      pass
    else:
      return match.group(1)

    try:
      # Different versions of raxml have different output :(
      (match,) = re.finditer("Final tree written to:\s+(\S+)",
                             raxml_stdout)
    except ValueError:
      # Didn't get one match again, assume something has gone wrong
      return None
    else:
      return match.group(1)

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
