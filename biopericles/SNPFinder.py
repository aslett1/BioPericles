import os
import re

from vcf import Reader
from biopericles.Common import LoadFastaMixin, \
                               RunExternalApplicationMixin, \
                               ExternalApplicationException

class SNPSitesException(ExternalApplicationException):
  pass

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

class SNPFinder(LoadFastaMixin, RunExternalApplicationMixin):
  def __init__(self):
    self.sequences = {}
    self.snps = None

  def add_snps(self):
    pass

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
