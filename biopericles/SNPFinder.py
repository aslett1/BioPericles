import Bio.SeqIO
import logging
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
    """Replace 'AB' info with 'GT' format

    If the FORMAT is '.' and the INFO is 'AB' change the FORMAT to be 'GT'
    otherwise don't change the line"""
    row = re.split(self._separator, line)
    ref = row[3]
    alt = row[4]
    fmt = row[8]
    info = row[7]
    gt_lookup = {base: str(index+1) for index,base in enumerate(alt.split(','))}
    gt_lookup['.'] = '0'
    if fmt == '.' and info == 'AB':
      row[7] = '.'
      row[8] = 'GT' # Set the format of this line to Genotype (GT)
      self._amend_genotype_fields(row, gt_lookup)
      return "\t".join(row)
    else:
      return line

  def _amend_genotype_fields(self, row, gt_lookup):
    new_fields = map(lambda field: gt_lookup.get(field, '.'), row[9:])
    row[9:] = new_fields

  def next(self):
    """Wraps PyVCF Reader's next method to ammend the relevant lines"""
    reader = self.reader
    self.reader = (self._amend_line(line) for line in reader)
    record = super(SNPSitesReader, self).next()
    self.reader = reader
    return record

class DeletableFile(object):
  def __init__(self, file_handle):
    """Wraps a file like object to signify that it's deletable

    We don't want to delete the user's files so we wrap ambiguous files in this
    object to signify that they were locally created and therefore can be
    removed"""
    self.file_handle = file_handle

  def remove(self):
    try:
      self.file_handle.close()
      os.remove(self.file_handle.name)
    except (OSError, AttributeError):
      # The file didn't exist or had never been created
      pass

  def __iter__(self):
    return self.file_handle

  def __getattr__(self, key):
    return getattr(self.file_handle, key)

class SNPFeatureBuilder(LoadFastaMixin, RunExternalApplicationMixin):
  def __init__(self):
    self.sequences = {}
    self.vcf_input_file = DeletableFile(None)
    self.feature_labels = [] # list of feature label objects
    self.features = {} # map of sample names to feature objects
    self.snp_sites_exec = "snp-sites"
    self.logger = logging.getLogger(__name__)

  def __del__(self):
    try:
      self.vcf_input_file.remove()
    except AttributeError:
      # The file might have already been deleted or might not have been created
      # by us so we shouldn't delete it.
      pass

  def create_vcf_from_sequences(self):
    with context_aware_tempfile('w', delete=False) as temp_fasta_file:
      try:
        self.vcf_input_file.remove()
      except AttributeError:
        # The file might have already been deleted or might not have been created
        # by us so we shouldn't delete it.
        pass
      self.vcf_input_file = DeletableFile(tempfile.NamedTemporaryFile('w', delete=False))
      self.logger.info("Writing sequences to %s" % temp_fasta_file.name)
      self._write_sequences(self.sequences.values(), temp_fasta_file)
      temp_fasta_file.close()
      self.vcf_input_file.close()

      with context_aware_tempdir() as output_directory:
        self.logger.info("Writing snps to %s" % self.vcf_input_file.name)
        self._run_snp_sites(self.snp_sites_exec, {'-o': self.vcf_input_file.name},
                            temp_fasta_file.name, output_directory)

      self.vcf_input_file = DeletableFile(open(self.vcf_input_file.name, 'r'))

  def create_features(self):
    self.features = {}
    self.feature_labels = []
    self.logger.info("Creating features")
    for record in self._get_records_from_vcf():
      self._add_record_to_features(record)

  def _get_records_from_vcf(self):
    self.vcf_input_file.seek(0)
    self.logger.info("Reading VCF records from %s" % self.vcf_input_file.name)
    return SNPSitesReader(self.vcf_input_file)

  def _add_record_to_features(self, record):
    feature_name = "SNP:{chromosome}:{position}".format(chromosome=record.CHROM,
                                                        position=record.POS)
    self.logger.debug("Adding feature '%s'" % feature_name)
    feature_updates = []
    try:
      for sample in record.samples:
        sample_name = sample.sample
        sample_value = 0 if sample.data.GT == '0' else 1
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
