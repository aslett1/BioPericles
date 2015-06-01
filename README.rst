BioPericles
===========

.. image:: https://travis-ci.org/sanger-pathogens/BioPericles.svg
    :target: https://travis-ci.org/sanger-pathogens/BioPericles

A work in progress. BioPericles is designed to use supervised Machine Learning
techniques to classify samples into clusters.  It currently extracts features 
using SNP analysis although other features could be used.

In the future it will hopefully be able to identfy which features are the
most descriminating between clusters.  This could be useful for creating
diagnostic tests.

Setup
-----

::

  python setup.py develop
  pip install -r requirements.txt

Scripts
-------

**cluster-sequences**

Takes a `multifasta <http://en.wikipedia.org/wiki/FASTA_format>`_ of sequences
and a metadata file mapping sequences to clusters and splits the file into one
file per cluster.

Requires: biopython

::

  export data=biopericles/tests/data
  cluster-sequences -o ${data} ${data}/clusters_spreadsheet.csv ${data}/multifasta.aln

**get-cluster-consensus**

Takes a list of multifastas, one per cluster.  For each multifasta it
calculates the bases which are consistent throughout the cluster and
marks those that are not with an 'N'.

Requires: biopython

::

  export data=biopericles/tests/data
  get-cluster-consensus -o ${data}/cluster-consensus.aln ${data}/multifasta.aln

**calculate-ancestral-sequences**

Takes an aligned multifasta, builds a phylogenetic tree from the sequences
(using `RAxML <https://github.com/stamatak/standard-RAxML>`_) and then
calculates what the sequences for the internal nodes might look like
(using `FASTML <http://fastml.tau.ac.il/>`_).  Outputs a fasta file and
newick tree to specified files or defaults to the current directory.

Requires: biopython, PyVCF, RAxML, FASTML

::

  export data=biopericles/tests/data
  calculate-ancestral-sequences ${data}/animals.mfa

**get-snp-features**

Takes an aligned multifasta and runs `snp-sites <https://github.com/sanger-pathogens/snp_sites>`_
to find `SNPs <http://en.wikipedia.org/wiki/Single-nucleotide_polymorphism>`_.
These SNPs are used to output a CSV with a row per sample and columns for each
SNP position.  Positions which match the reference (the first sequence in the
multifasta) have the value :code:`0` while differences are marked with :code:`1`.

You can also pass a `VCF file <http://en.wikipedia.org/wiki/Variant_Call_Format>`_
as input by setting the :code:`--input_format` option to
:code:`vcf`.  This bipasses any call to snp-sites.  VCFs can have a couple of
formats.  The first is that output by snp-sites which has :code:`INFO` set to
:code:`AB`, genotype fields set to either :code:`.` if this sample matches the
'reference' or the alternative base otherwise.  See
`file_with_SNPs.aln.vcf <biopericles/tests/data/file_with_SNPs.aln.vcf>`_ for an
example.  Alternatively you can provide a VCF with a Genotype field like
`file_with_SNPs_in_GT_format.aln.vcf <biopericles/tests/data/file_with_SNPs_in_GT_format.aln.vcf>`_.

For convenience, if you provide the script with a multifasta, you can ask to save
the SNPs calculated by snp-sites in VCF format using the :code:`--vcf-output` option.
This file is output in the Genotype style rather than as produced by snp-sites.
Users can then filter / edit this file before recalculating different features.

Requires: biopython, snp-sites, PyVCF

::

  export data=biopericles/tests/data
  get-snp-features ${data}/file_with_SNPs.aln
  # or alternatively
  get-snp-features --input_format vcf ${data}/file_with_SNPs_in_GT_format.aln.vcf

**train-classifier**

Takes a CSV of features and a CSV with sample to label mappings and trains a
classifier.  The features should be supplied with one row per sample and a
column per feature with :code:`1` or :code:`0` in each cell.

Although designed to work with the output from `get-snp-features <scripts/get-snp-features>`_
there is no reason you cannot add or remove features by altering the columns.

It is possible to pass in a second feature file using the :code:`--cluster-features`
option.  This file is used to get a list of features which should be used for
classifier training.  This is useful if you only want to use SNPs which are consistent
within clusters.

The classifier is output as a `Python Pickle <https://docs.python.org/2/library/pickle.html>`_
although this could be changed.  Pickles of unknown origin should be considered
dangerous; the load method on `SampleClassifier <biopericles/SampleClassifier.py>`_
expects to be passed an MD5 or SHA1 before loading a classifier from disc.  These
values are output at the end of classifier training.

Requires: numpy, scikit-learn, scipy

::

  export data=biopericles/tests/data
  train-classifier ${data}/fruit_example_features.csv ${data}/fruit_example_metadata.csv fruit_classifier.pkl

Prototypes
----------

The `scripts <scripts/>`_ are pretty reusable, `prototypes <prototype/>`_ are less polished
but might be a useful starting point for someone trying to use BioPericles on their own
data.  They're not very well tested so use them with care.

**annotate_vcf.py**

Takes a VCF, `GFF <http://en.wikipedia.org/wiki/General_feature_format>`_ and a Fasta file.
It uses these to output a VCF with anotations taken from the GFF file to identify
synonymous, nonsynonymous and intergenic variants by looking for `CDS <http://en.wikipedia.org/wiki/Coding_region>`_
features.  It uses a format similar to `Variant Effect Predictor <http://www.ensembl.org/info/docs/tools/vep/index.html>`_
but provides less detail.

Requires: `GenomeTools <http://genometools.org/>`_

**install_dependencies.sh**

Installing the dependencies for BioPericles can be a bit of a pain.  A script is
included in the prototypes directory to help install dependencies in a clean system
on which you have root access.  NB This script manipulates your PATH and 
LD_LIBRARY_PATH environment variables.  It has only been tested in an Ubuntu
Docker Container.

`install_dependencies.sh <prototype/install_dependencies.sh>`_ works a bit like a
Makefile.  As long as you run it from the same location it will not re-download or
build the dependencies and will just re-update your environment variables if required.
This is useful for Travis CI because it allows us to cache the build folder to speed up
subsequent test runs.

::

  cd BioPericles/prototype
  sudo ./install_dependencies.sh
  sudo pip install -r ../requirements.txt

TODO
----

- Use the classifier to classify things
- Make it clearer what the test set is or possible to specify what it should be
- Use the classifier to identify a subset of features which are the most discriminative
- Consider how the classifier could be tuned (either manually or automatically)
- Make it easier to score the classifier using a separate test set
