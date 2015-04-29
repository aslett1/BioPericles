BioPericles
===========

A work in progress.  At the moment it can split a multifasta file using cluster data in a CSV.

Examples
--------

**Setup**

::

  python setup.py develop


**cluster_sequences**

::

  export data=biopericles/tests/data
  cluster-sequences -o ${data} ${data}/clusters_spreadsheet.csv ${data}/multifasta.aln

**get-cluster-consensus**

Takes a list of multifastas, one per cluster.  For each multifasta it calculates
the bases which are consistent throughout the cluster and marks those that are 
not with an 'N'.

NB This script currently requires at least ten times as much memory as the size
of the biggest input file passed to it.

::

  export data=biopericles/tests/data
  get-cluster-consensus -o ${data}/cluster-consensus.aln ${data}/multifasta.aln
