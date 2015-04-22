import os
import shutil
import tempfile
import unittest

from biopericles.SNPFinder import SNPFinder

def test_data():
  this_file = os.path.abspath(__file__)
  this_dir = os.path.dirname(this_file)
  return os.path.join(this_dir, 'data')

class TestSNPFinder(unittest.TestCase):
  def test_run_snp_sites(self):
    finder = SNPFinder()

    fasta_filename = os.path.join(test_data(), 'file_with_SNPs.aln')
    output_directory = tempfile.mkdtemp()

    stdout, stderr = finder._run_snp_sites('snp-sites', {},
                                           fasta_filename,
                                           output_directory)

    expected_output_filename = os.path.join(output_directory, 'all_snps.vcf')
    self.assertTrue(os.path.isfile(expected_output_filename))

    shutil.rmtree(output_directory)
