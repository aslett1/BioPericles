import Bio.SeqIO
import os
import re
import shutil
import subprocess
import tempfile

from collections import OrderedDict
from copy import copy

def try_and_get_filename(filehandle):
  try:
    return filehandle.name
  except AttributeError:
    return "ANONYMOUS_FILE"

class RaxmlException(Exception):
  def __init__(self, message, returncode, stdout, stderr):
    super(RaxmlException, self).__init__(message)
    self.returncode = returncode
    self.stdout = stdout
    self.stderr = stderr

class FastmlException(Exception):
  def __init__(self, message, returncode, stdout, stderr):
    super(FastmlException, self).__init__(message)
    self.returncode = returncode
    self.stdout = stdout
    self.stderr = stderr

class TreeBuilder(object):
  def __init__(self):
    self.sequences = None # a dictionary of {sequence_name: SeqIO object}
    self.tree = None # a BioPython Phylo tree

  def load_fasta_sequences(self, fasta_file):
    """Load sequences from a fasta_file into a dictionary of {sequence_name: SeqIO object}

    Takes a filehandle to a aligned multifasta
    """

    list_of_sequences = Bio.SeqIO.parse(fasta_file, 'fasta')
    self.sequences = OrderedDict()
    for seq in list_of_sequences:
      if seq.name in self.sequences:
        message = "Tried to load sequences from {filename} but got multiple sequences for {sequence}"
        raise ValueError(message.format(filename=try_and_get_filename(fasta_file),
                                        sequence=seq.name))
      self.sequences[seq.name] = seq

  def build_tree(self):
    phylip = self._create_temporary_phylip(self.sequences)
    phylip.close()
    output_directory = tempfile.mkdtemp()
    raxml_stdout, raxml_stderr = self._run_raxml('raxmlHPC', {}, phylip.name,
                                                 output_directory)
    raxml_tree_filename = self._get_raxml_tree_file(raxml_stdout)
    (tree,) = Bio.Phylo.parse(raxml_tree_filename, 'newick')
    tree.root_at_midpoint()
    os.remove(phylip.name)
    shutil.rmtree(output_directory)
    return tree

  def add_hereditary_nodes(self, tree, sequences):
    pass # takes a tree, a dictionary of BioPython sequences; returns a tree object

  def write_output(self, output_directory, filename_prefix):
    pass # writes the sequences and internal nodes and full tree

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
    arguments_dict = self._merge_commandline_arguments(default_arguments,
                                                       raxml_arguments)
    arguments_list = self._build_commandline_arguments(arguments_dict)
    raxml_process = subprocess.Popen([raxml_executable] + arguments_list,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE
                                    )
    raxml_stdout, raxml_stderr = raxml_process.communicate()
    if raxml_process.returncode != 0:
      raise RaxmlException("Problem running raxml on %s; some output in %s" %
                         (phylip_filename, output_directory),
                           raxml_process.returncode,
                           raxml_stdout,
                           raxml_stderr)
    return (raxml_stdout, raxml_stderr)

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

    arguments_dict = self._merge_commandline_arguments(default_arguments,
                                                       fastml_arguments)
    arguments_list = self._build_commandline_arguments(arguments_dict)
    fastml_process = subprocess.Popen([fastml_executable] + arguments_list,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     cwd=output_directory
                                    )
    fastml_stdout, fastml_stderr = fastml_process.communicate()
    if fastml_process.returncode != 0:
      raise FastmlException("Problem running fastml using %s and %s; some output in %s" %
                           (sequence_fasta_filename, tree_filename, output_directory),
                           fastml_process.returncode,
                           fastml_stdout,
                           fastml_stderr)
    return (fastml_stdout, fastml_stderr)

  def _merge_commandline_arguments(self, default_arguments, new_arguments):
    """Takes a dictionary mapping strings to strings and merges it with another

    If the value is None, it removes that key, value from the output;
    If the key is in the default_arguments it is overwriden by new_arguments;
    If the key isn't in default_arguments, it is added"""
    arguments = copy(default_arguments)
    for key,value in new_arguments.items():
      if value == None and key in arguments:
        del arguments[key]
      elif value != None:
        arguments[key] = value
    return arguments

  def _build_commandline_arguments(self, arguments):
    """Turns an ordered dict of arguments into a list

    Ignores values which are '';
    Ignores keys and values where the value is None"""
    output = []
    for key,value in arguments.items():
      if value == None:
        continue
      elif value == '':
        output.append(key)
      else:
        output += [key, value]

    return output
