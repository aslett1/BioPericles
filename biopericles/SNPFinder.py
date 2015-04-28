import Bio.SeqIO
import os
import re
import shutil
import tempfile

from vcf import Reader
from biopericles.Common import LoadFastaMixin, \
                               RunExternalApplicationMixin, \
                               ExternalApplicationException

class SNPSitesException(ExternalApplicationException):
  pass

class SNPFeature(object):
  def __init__(self, name, label, features):
    self.name = name
    self.label = label
    self.features = features

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
    self.vcf = None
    self.vcf_output_file = None
    self.temp_vcf_file = None

  def add_vcf(self):
    temp_fasta_file = tempfile.NamedTemporaryFile('w', delete=False)
    self.temp_vcf_file = tempfile.NamedTemporaryFile('w', delete=False)
    self._write_sequences(self.sequences.values(), temp_fasta_file)
    temp_fasta_file.close()
    self.temp_vcf_file.close()

    output_directory = tempfile.mkdtemp()

    self._run_snp_sites('snp-sites', {'-o': self.temp_vcf_file.name},
                        temp_fasta_file.name, output_directory)

    self.temp_vcf_file = open(self.temp_vcf_file.name, 'r')
    self.vcf = SNPSitesReader(self.temp_vcf_file)

    os.remove(temp_fasta_file.name)
    shutil.rmtree(output_directory)

  def create_features(self):
    pass

  def write_vcf(self):
    pass

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
