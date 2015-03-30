BioPericles
===========

A work in progress.  At the moment it can split a multifasta file using cluster data in a CSV.

::

  python setup.py develop
  export data=biopericles/tests/data
  cluster-sequences -o ${data} ${data}/clusters_spreadsheet.csv ${data}/multifasta.aln
