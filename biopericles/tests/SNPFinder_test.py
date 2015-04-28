import os
import shutil
import tempfile
import unittest
import vcf

from mock import MagicMock
from StringIO import StringIO
from biopericles.SNPFinder import SNPFeatureBuilder, SNPSitesReader

def test_data():
  this_file = os.path.abspath(__file__)
  this_dir = os.path.dirname(this_file)
  return os.path.join(this_dir, 'data')

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
  def test_add_vcf(self):
    pass

  def test_create_features(self):
    pass

  def test_write_vcf(self):
    pass

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
