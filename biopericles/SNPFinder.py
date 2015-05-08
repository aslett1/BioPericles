import Bio.SeqIO
import os
import re
import shutil
import tempfile

from vcf import Reader
from biopericles.Common import LoadFastaMixin, \
                               RunExternalApplicationMixin, \
                               ExternalApplicationException, \
                               context_aware_tempfile, context_aware_tempdir

class SNPSitesException(ExternalApplicationException):
  pass

class SNPFeature(object):
  def __init__(self, name, label, features):
    self.name = name
    self.label = label
    self.features = features

class SNPFeatureLabel(object):
  def __init__(self, chromosome, position):
    self.chromosome = chromosome
    self.position = position

class SNPSitesReader(Reader):
  """Wraps PyVCF to make it compatible with our VCF Files

  The VCF files we generally use will have a format of '.' and simply provide
  '.' for each sample if it shares a common base or the alternative base if it
  is a SNP.  PyVCF doesn't like this format so we monkey patch the file reader
  to ammend the relevant lines before PyVCF has a chance to parse them."""
  def _amend_line(self, line):
    """Fix the format for this line

    If the FORMAT is '.' and the INFO is 'AB' change the FORMAT to also be
    'AB'"""
    row = re.split(self._separator, line)
    fmt = row[8]
    info = row[7]
    if fmt == '.' and info == 'AB':
      row[8] = 'AB' # Set the format of this line to Alternative Base (AB)
      return "\t".join(row)
    else:
      return line

  def next(self):
    """Wraps PyVCF Reader's next method to ammend the relevant lines"""
    reader = self.reader
    self.reader = (self._amend_line(line) for line in reader)
    record = super(SNPSitesReader, self).next()
    self.reader = reader
    return record

class SNPFeatureBuilder(LoadFastaMixin, RunExternalApplicationMixin):
  def __init__(self):
    self.sequences = {}
    self.vcf_output_file = None
    self.vcf_input_file = None
    self.feature_labels = [] # list of feature label objects
    self.features = {} # map of sample names to feature objects
    self.snp_sites_exec = "snp-sites"

  def __del__(self):
    try:
      self.vcf_input_file.close()
      os.remove(self.vcf_input_file.name)
    except (OSError, AttributeError):
      # The file didn't exist or had never been created
      pass

  def create_vcf_from_sequences(self):
    with context_aware_tempfile('w', delete=False) as temp_fasta_file:
      if self.vcf_input_file != None:
        self.vcf_input_file.close()
        os.remove(self.vcf_input_file.name)
      self.vcf_input_file = tempfile.NamedTemporaryFile('w', delete=False)
      self._write_sequences(self.sequences.values(), temp_fasta_file)
      temp_fasta_file.close()
      self.vcf_input_file.close()

      with context_aware_tempdir() as output_directory:
        self._run_snp_sites(self.snp_sites_exec, {'-o': self.vcf_input_file.name},
                            temp_fasta_file.name, output_directory)

      self.vcf_input_file = open(self.vcf_input_file.name, 'r')

  def create_features(self):
    self.features = {}
    self.feature_labels = []
    for record in self._get_records_from_vcf():
      self._add_record_to_features(record)

  def _get_records_from_vcf(self):
    self.vcf_input_file.seek(0)
    return SNPSitesReader(self.vcf_input_file)

  def _add_record_to_features(self, record):
    feature_name = "SNP:{chromosome}:{position}".format(chromosome=record.CHROM,
                                                        position=record.POS)
    feature_updates = []
    try:
      for sample in record.samples:
        sample_name = sample.sample
        sample_value = 0 if sample.data.AB == None else 1
        feature_updates.append((sample_name, sample_value))
    except AttributeError:
      # One or more of the samples didn't have alternative base ('AB') data so
      # we should ignore this record without updating it.
      return

    for sample_name, sample_value in feature_updates:
      self.features.setdefault(sample_name, []).append(sample_value)

    self.feature_labels.append(feature_name)

  def _write_sequences(self, sequences, output_file):
    Bio.SeqIO.write(sequences, output_file, 'fasta')

  def _run_snp_sites(self, snp_sites_executable, arguments, fasta_filename,
                     output_directory):
    default_arguments = {
      "-v": "",
      "-o": os.path.join(output_directory, 'all_snps.vcf'),
      fasta_filename: ""
    }
    stdout, stderr, returncode = self._run_application(snp_sites_executable,
                                                       default_arguments,
                                                       arguments,
                                                       cwd=output_directory)
    if returncode != 0:
      raise SNPSitesException("Problem running snp_sites on %s; some output in %s" %
                              (fasta_filename, output_directory),
                              returncode,
                              stdout,
                              stderr)
    return (stdout, stderr)
