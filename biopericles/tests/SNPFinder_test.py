import os
import shutil
import tempfile
import unittest
import vcf

from mock import MagicMock, patch
from StringIO import StringIO

from biopericles.SNPFinder import SNPFeatureBuilder, SNPSitesReader

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

    line = "0\t1\t2\t3\t4\t5\t6\tAB\t.\t9\t10"
    expected = "0\t1\t2\t3\t4\t5\t6\tAB\tAB\t9\t10"
    self.assertEqual(reader._amend_line(line), expected)

    # It would be elegant if it maintained the separator but this is unlikely to
    # case big issues, hopefully
    line = "0 1 2 3 4 5 6 AB . 9 10"
    expected = "0\t1\t2\t3\t4\t5\t6\tAB\tAB\t9\t10"
    self.assertEqual(reader._amend_line(line), expected)

    line = "0  1 2    3 4\t5\t6 AB . 9 10"
    expected = "0\t1\t2\t3\t4\t5\t6\tAB\tAB\t9\t10"
    self.assertEqual(reader._amend_line(line), expected)

    line = "0\t1\t2\t3\t4\t5\t6\tAB\tGT\t9\t10"
    expected = "0\t1\t2\t3\t4\t5\t6\tAB\tGT\t9\t10"
    self.assertEqual(reader._amend_line(line), expected)
    SNPSitesReader.__bases__ = snp_sites_old_bases

  def test_next(self):
    vcf_filename = os.path.join(test_data(), 'file_with_SNPs.aln.vcf')
    vcf_file = open(vcf_filename, 'r')

    reader = SNPSitesReader(vcf_file)
    record = reader.next()

    samples = record.samples
    samples_with_alternative_bases = [sample.sample for sample in samples if
                                      sample.data.AB]

    expected = ['3002_8_1', '3002_8_2', '3002_8_6', '4056_2_10', '4056_2_4',
                '4056_8_6', '5174_5_1', '5174_5_7', '5174_5_9', '5174_6_10',
                '5174_7_1', '5174_8_5']

    self.assertItemsEqual(samples_with_alternative_bases, expected)

class TestSNPFeatureBuilder(unittest.TestCase):
  @patch("biopericles.SNPFinder.tempfile")
  def test_create_vcf_from_sequences(self, temp_mock):
    builder = SNPFeatureBuilder()

    temp_fasta_file = tempfile.NamedTemporaryFile('w', delete=False)
    temp_vcf_file = tempfile.NamedTemporaryFile('w', delete=False)
    random_folder = tempfile.mkdtemp()

    temp_mock.NamedTemporaryFile.side_effect = [temp_fasta_file, temp_vcf_file]
    temp_mock.mkdtemp.return_value = random_folder

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    fasta_file = open(fasta_filename, 'r')

    builder.load_fasta_sequences(fasta_file)
    builder.create_vcf_from_sequences()

    self.assertFalse(os.path.isfile(temp_fasta_file.name))
    self.assertFalse(os.path.isdir(random_folder))
    self.assertTrue(os.path.isfile(temp_vcf_file.name))

    builder.vcf_input_file.seek(0)
    records = SNPSitesReader(builder.vcf_input_file)
    number_of_records = sum((1 for record in records))

    self.assertEqual(number_of_records, 5)

    fasta_file.close()
    temp_vcf_file.close()
    os.remove(temp_vcf_file.name)

  @patch("biopericles.SNPFinder.os.remove")
  def test_create_vcf_from_sequences_deletes_existing(self, remove_mock):
    builder = SNPFeatureBuilder()
    temp_vcf_file = tempfile.NamedTemporaryFile('w', delete=True)
    builder.vcf_input_file = temp_vcf_file

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    fasta_file = open(fasta_filename, 'r')

    builder.load_fasta_sequences(fasta_file)
    builder.create_vcf_from_sequences()

    remove_mock.assert_any_call(temp_vcf_file.name)

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
    output_directory = tempfile.mkdtemp()

    stdout, stderr = builder._run_snp_sites('snp-sites', {},
                                           fasta_filename,
                                           output_directory)

    expected_output_filename = os.path.join(output_directory, 'all_snps.vcf')
    self.assertTrue(os.path.isfile(expected_output_filename))

    shutil.rmtree(output_directory)
