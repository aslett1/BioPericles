import Bio.SeqIO
import logging
import subprocess

from collections import OrderedDict
from copy import copy

def try_and_get_filename(filehandle):
  try:
    return filehandle.name
  except AttributeError:
    return "ANONYMOUS_FILE"

class ExternalApplicationException(Exception):
  def __init__(self, message, returncode, stdout, stderr):
    super(ExternalApplicationException, self).__init__(message)
    self.returncode = returncode
    self.stdout = stdout
    self.stderr = stderr

class LoadFastaMixin(object):
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
    return self.sequences

class RunExternalApplicationMixin(object):
  def _run_application(self, executable, default_arguments,
                       additional_arguments, **kwargs):
    logger = logging.getLogger(__name__)
    arguments_dict = self._merge_commandline_arguments(default_arguments,
                                                       additional_arguments)
    arguments_list = self._build_commandline_arguments(arguments_dict)
    process = subprocess.Popen([executable] + arguments_list,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               **kwargs
                              )
    logger.debug("Called external application '%s'" % " ".join([executable] +
                                                               arguments_list))
    stdout, stderr = process.communicate()
    return (stdout, stderr, process.returncode)

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

