import os

from biopericles.Common import LoadFastaMixin, \
                               RunExternalApplicationMixin, \
                               ExternalApplicationException

class SNPSitesException(ExternalApplicationException):
  pass

class SNPFinder(LoadFastaMixin, RunExternalApplicationMixin):
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
