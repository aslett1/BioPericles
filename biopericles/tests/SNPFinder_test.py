import os
import shutil
import tempfile
import unittest
import vcf

from contextlib import contextmanager
from mock import MagicMock, patch
from StringIO import StringIO

from biopericles.SNPFinder import SNPFeatureBuilder, SNPSitesReader
from biopericles.Common import context_aware_tempfile, context_aware_tempdir

def test_data():
  this_file = os.path.abspath(__file__)
  this_dir = os.path.dirname(this_file)
  return os.path.join(this_dir, 'data')

class FakeRecord(object):
  """Creates 'Records' similar to those parsed from vcf files

  Each vcf file is normally parsed into records, one per SNP.  This object is
  similar enough for tests"""
  def __init__(self, chromosome, position, calls):
    self.CHROM = chromosome
    self.POS = position
    sample_names = ["sample_%s" % i for i in range(len(calls))]
    self.samples = [FakeCall(*args) for args in zip(sample_names, calls)]

class FakeCall(object):
  """Records have a list of 'Calls' for each sample"""
  def __init__(self, sample, call):
    self.data = FakeDatum(call)
    self.sample = sample

class FakeDatum(object):
  """The calls have values, in this case we're only mocking the AB value"""
  def __init__(self, datum):
    self.AB = datum

class TestSNPSitesReader(unittest.TestCase):
  def test_ammend_line(self):
    snp_sites_old_bases = SNPSitesReader.__bases__
    SNPSitesReader.__bases__ = (MagicMock,) # I don't want to test vcf.Reader
    reader = SNPSitesReader()
    reader._separator="\t| +"

    line = "0\t1\t2\tC\tT,G\t5\t6\tAB\t.\tG\t.\tT"
    expected = "0\t1\t2\tC\tT,G\t5\t6\t.\tGT\t2\t0\t1"
    self.assertEqual(reader._amend_line(line), expected)

    # It would be elegant if it maintained the separator but this is unlikely to
    # case big issues, hopefully
    line = "0 1 2 C T,G 5 6 AB . G . T"
    expected = "0\t1\t2\tC\tT,G\t5\t6\t.\tGT\t2\t0\t1"
    self.assertEqual(reader._amend_line(line), expected)

    line = "0  1 2    C T,G\t5\t6 AB . G . T"
    expected = "0\t1\t2\tC\tT,G\t5\t6\t.\tGT\t2\t0\t1"
    self.assertEqual(reader._amend_line(line), expected)

    line = "0\t1\t2\tC\tT,G\t5\t6\tAB\tGT\t2\t0\t1"
    expected = "0\t1\t2\tC\tT,G\t5\t6\tAB\tGT\t2\t0\t1"
    self.assertEqual(reader._amend_line(line), expected)
    SNPSitesReader.__bases__ = snp_sites_old_bases

  def test_next(self):
    vcf_filename = os.path.join(test_data(), 'file_with_SNPs.aln.vcf')
    vcf_file = open(vcf_filename, 'r')

    reader = SNPSitesReader(vcf_file)
    record = reader.next()

    samples = record.samples
    samples_with_alternative_bases = [sample.sample for sample in samples if
                                      sample.data.GT != '0']

    expected = ['3002_8_1', '3002_8_2', '3002_8_6', '4056_2_10', '4056_2_4',
                '4056_8_6', '5174_5_1', '5174_5_7', '5174_5_9', '5174_6_10',
                '5174_7_1', '5174_8_5']

    self.assertItemsEqual(samples_with_alternative_bases, expected)

def create_context_aware_tempfile_mock(filenames):
  """Creates something like a context_aware_tempfile but remembers created files

  context_aware_tempfiles are deleted at the end of the context.  This does
  exactly the same thing but also populates the filenames array with details of
  the files created so that you can check they have been deleted"""
  @contextmanager
  def context_aware_tempfile_mock(*args, **kwargs):
    with context_aware_tempfile(*args, **kwargs) as tempfile:
      filenames.append(tempfile.name)
      yield tempfile
  return context_aware_tempfile_mock

def create_context_aware_tempdir_mock(folder_names):
  """Creates something like a context_aware_tempfile but remembers created files

  context_aware_tempfiles are deleted at the end of the context.  This does
  exactly the same thing but also populates the filenames array with details of
  the files created so that you can check they have been deleted"""
  @contextmanager
  def context_aware_tempdir_mock(*args, **kwargs):
    with context_aware_tempdir(*args, **kwargs) as tempdir:
      folder_names.append(tempdir)
      yield tempdir
  return context_aware_tempdir_mock

class TestSNPFeatureBuilder(unittest.TestCase):
  @patch("biopericles.SNPFinder.context_aware_tempdir")
  @patch("biopericles.SNPFinder.context_aware_tempfile")
  @patch("biopericles.SNPFinder.tempfile")
  def test_create_vcf_from_sequences(self, temp_mock, ctx_tempfile_mock,
                                     ctx_tempdir_mock):
    builder = SNPFeatureBuilder()

    temp_vcf_file = tempfile.NamedTemporaryFile('w', delete=False)

    temp_files = []
    temp_folders = []
    ctx_tempfile_mock.side_effect = create_context_aware_tempfile_mock(temp_files)
    ctx_tempdir_mock.side_effect = create_context_aware_tempdir_mock(temp_folders)
    temp_mock.NamedTemporaryFile.return_value = temp_vcf_file

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    fasta_file = open(fasta_filename, 'r')

    builder.load_fasta_sequences(fasta_file)
    builder.create_vcf_from_sequences()

    self.assertEqual(len(temp_files), 1)
    self.assertFalse(os.path.isfile(temp_files[0]))
    self.assertEqual(len(temp_folders), 1)
    self.assertFalse(os.path.isdir(temp_folders[0]))
    self.assertTrue(os.path.isfile(temp_vcf_file.name))

    builder.vcf_input_file.seek(0)
    records = SNPSitesReader(builder.vcf_input_file)
    number_of_records = sum((1 for record in records))

    self.assertEqual(number_of_records, 5)

    fasta_file.close()
    temp_vcf_file.close()
    os.remove(temp_vcf_file.name)

  @patch("biopericles.SNPFinder.os")
  def test_create_vcf_from_sequences_deletes_existing(self, os_mock):
    builder = SNPFeatureBuilder()
    os_mock.remove.side_effect = os.remove
    with context_aware_tempfile('w', delete=False) as temp_vcf_file:
      builder.vcf_input_file = temp_vcf_file

      fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
      fasta_file = open(fasta_filename, 'r')

      builder.load_fasta_sequences(fasta_file)
      builder.create_vcf_from_sequences()

      os_mock.remove.assert_any_call(temp_vcf_file.name)

  @patch("biopericles.SNPFinder.os")
  def test_create_vcf_from_sequences_deletes_when_done(self, os_mock):
    builder = SNPFeatureBuilder()
    os_mock.remove.side_effect = os.remove

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    fasta_file = open(fasta_filename, 'r')

    builder.load_fasta_sequences(fasta_file)
    builder.create_vcf_from_sequences()

    temp_vcf_filename = builder.vcf_input_file.name
    del(builder)
    os_mock.remove.assert_any_call(temp_vcf_filename)
    try:
      os.remove(temp_vcf_filename)
    except OSError:
      # Looks like it's already deleted
      pass

  def test_get_records_from_vcf(self):
    builder = SNPFeatureBuilder()

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    fasta_file = open(fasta_filename, 'r')

    builder.load_fasta_sequences(fasta_file)
    builder.create_vcf_from_sequences()

    records = builder._get_records_from_vcf()
    number_of_records = sum((1 for record in records))
    self.assertEqual(number_of_records, 5)

    records = builder._get_records_from_vcf()
    number_of_records = sum((1 for record in records))
    self.assertEqual(number_of_records, 5)

  def test_create_features(self):
    builder = SNPFeatureBuilder()
    fake_vcf = []
    fake_vcf.append(FakeRecord(1,1,[None, None, 'A']))
    fake_vcf.append(FakeRecord(1,2,[None, 'G', 'T']))
    fake_vcf.append(FakeRecord(1,3,['A', 'C', None]))
    builder._get_records_from_vcf=MagicMock(return_value=fake_vcf)

    builder.create_features()

    self.assertItemsEqual(builder.features.keys(), ['sample_0', 'sample_1',
                                                    'sample_2'])
    self.assertEqual(builder.features['sample_0'], [0,0,1])
    self.assertEqual(builder.features['sample_1'], [0,1,1])
    self.assertEqual(builder.features['sample_2'], [1,1,0])
    self.assertEqual(builder.feature_labels, ['SNP:1:1', 'SNP:1:2', 'SNP:1:3'])

    fake_vcf = [FakeRecord(2,1,['T', None])]
    builder._get_records_from_vcf=MagicMock(return_value=fake_vcf)

    builder.create_features()

    self.assertItemsEqual(builder.features.keys(), ['sample_0', 'sample_1'])
    self.assertEqual(builder.features['sample_0'], [1])
    self.assertEqual(builder.features['sample_1'], [0])
    self.assertEqual(builder.feature_labels, ['SNP:2:1'])

  def test_add_record_to_features(self):
    builder = SNPFeatureBuilder()

    record = FakeRecord(1,1,[None, None, 'A'])
    builder._add_record_to_features(record)

    self.assertItemsEqual(builder.features.keys(), ['sample_0', 'sample_1',
                                                    'sample_2'])
    self.assertEqual(builder.features['sample_0'], [0])
    self.assertEqual(builder.features['sample_1'], [0])
    self.assertEqual(builder.features['sample_2'], [1])
    self.assertEqual(builder.feature_labels, ['SNP:1:1'])

    record = FakeRecord(1,2,[None, 'G', 'T'])
    builder._add_record_to_features(record)

    self.assertItemsEqual(builder.features.keys(), ['sample_0', 'sample_1',
                                                    'sample_2'])
    self.assertEqual(builder.features['sample_0'], [0,0])
    self.assertEqual(builder.features['sample_1'], [0,1])
    self.assertEqual(builder.features['sample_2'], [1,1])
    self.assertEqual(builder.feature_labels, ['SNP:1:1', 'SNP:1:2'])

    record = FakeRecord(1,3,['A', 'C', None])
    record.samples[1].data.__delattr__('AB') # One of the samples is missing an 'AB'
    builder._add_record_to_features(record)

    self.assertItemsEqual(builder.features.keys(), ['sample_0', 'sample_1',
                                                    'sample_2'])
    self.assertEqual(builder.features['sample_0'], [0,0])
    self.assertEqual(builder.features['sample_1'], [0,1])
    self.assertEqual(builder.features['sample_2'], [1,1])
    self.assertEqual(builder.feature_labels, ['SNP:1:1', 'SNP:1:2'])

  def test_run_snp_sites(self):
    builder = SNPFeatureBuilder()

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')

    with context_aware_tempdir() as output_directory:
      stdout, stderr = builder._run_snp_sites('snp-sites', {},
                                             fasta_filename,
                                             output_directory)

      expected_output_filename = os.path.join(output_directory, 'all_snps.vcf')
      self.assertTrue(os.path.isfile(expected_output_filename))
