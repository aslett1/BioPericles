import os
import random
import re
import shutil
import tempfile
import unittest
import Bio.Phylo

from mock import patch, MagicMock
from StringIO import StringIO

from biopericles.TreeBuilder import TreeBuilder, RaxmlException, FastmlException

def test_data():
  this_file = os.path.abspath(__file__)
  this_dir = os.path.dirname(this_file)
  return os.path.join(this_dir, 'data')

def create_fasta_without_comments(filename):
  """Create a fasta file without comments

  This is probably not the most efficient way but at least it is cross platform
  and it runs plenty quickly for these tests."""
  fasta_file_without_comments = tempfile.NamedTemporaryFile(delete=False)
  with open(filename, 'r') as original_fasta_file:
    for line in original_fasta_file:
      line_without_comments = re.sub(r'\s*;.*$', '', line)
      if line_without_comments.strip() != '':
        fasta_file_without_comments.write(line_without_comments)
  return fasta_file_without_comments.name

class TestTreeBuilder(unittest.TestCase):
  @patch('biopericles.TreeBuilder.tempfile.NamedTemporaryFile')
  def test_create_temporary_phylip(self, temp_mock):
    temp_file = StringIO()
    temp_mock.return_value = temp_file

    builder = TreeBuilder()

    fasta_file = StringIO("""\
>cluster_A
NNAACAAAANN
>cluster_B
CCAACAAAANN
""")

    expected_output = """\
 2 11
cluster_A  NNAACAAAAN N
cluster_B  CCAACAAAAN N
"""

    builder.load_fasta_sequences(fasta_file)
    output_file = builder._create_temporary_phylip(builder.sequences)

    self.assertEqual(output_file, temp_file)

    output_file.seek(0)
    self.assertEqual(output_file.read(), expected_output)

  @patch("biopericles.TreeBuilder.tempfile")
  def test_build_tree(self, temp_mock):
    builder = TreeBuilder()

    # I could mock _run_raxml but, because I'm using a regex to parse the
    # results, I want to know straight away if the output of RAxML changes in a
    # way that breaks this.

    random_file = tempfile.NamedTemporaryFile('w', delete=False)
    random_filename = random_file.name
    random_folder = tempfile.mkdtemp()

    temp_mock.NamedTemporaryFile.return_value = random_file
    temp_mock.mkdtemp.return_value = random_folder

    fasta_file = open(os.path.join(test_data(), 'animals.mfa'), 'r')
    builder.load_fasta_sequences(fasta_file)
    fasta_file.close()

    tree = builder.build_tree()

    tree_nodes = sorted([node.name for node in tree.get_terminals()])
    expected_tree_nodes = ["Carp", "Chicken", "Cow", "Frog", "Human", "Loach",
                           "Mouse", "Rat", "Seal", "Whale"]

    self.assertIsInstance(tree, Bio.Phylo.Newick.Tree)
    self.assertEqual(tree_nodes, expected_tree_nodes)
    self.assertFalse(os.path.isdir(random_folder))
    self.assertFalse(os.path.isfile(random_filename))

  def test_run_raxml(self):
    builder = TreeBuilder()

    fasta_file = open(os.path.join(test_data(), 'animals.mfa'), 'r')
    builder.load_fasta_sequences(fasta_file)
    fasta_file.close()

    phylip_file = builder._create_temporary_phylip(builder.sequences)
    phylip_file.close()

    output_directory = tempfile.mkdtemp()
    raxml_stdout, raxml_stderr = builder._run_raxml('raxmlHPC', {}, phylip_file.name, output_directory)

    self.assertTrue('Best-scoring ML tree written to' in raxml_stdout)

    os.remove(phylip_file.name)
    shutil.rmtree(output_directory)

  def test_get_raxml_tree_file(self):
    builder = TreeBuilder()

    fake_stdout = """\
Inference[0] final GAMMA-based Likelihood: -385.279958 tree written to file/tmp/tmpaA9uJN/RAxML_result.raxml_output


Starting final GAMMA-basedthorough Optimization on tree 0 likelihood -385.279958 .... 

FinalGAMMA-based Score of best tree -385.279958

Program execution info written to/tmp/tmpaA9uJN/RAxML_info.raxml_output
Best-scoring ML tree written to:  /tmp/tmpaA9uJN/RAxML_bestTree.raxml_output

Overall execution time: 0.219131  secs or 0.000061 hours or 0.000003 days

"""

    output = builder._get_raxml_tree_file(fake_stdout)
    self.assertEqual(output, '/tmp/tmpaA9uJN/RAxML_bestTree.raxml_output')

    fake_stdout = """\
Oh dear

Something went wrong

I've written something to /tmp/wrong.log

But it is definitly not what you're looking for.

Sorry
"""

    output = builder._get_raxml_tree_file(fake_stdout)
    self.assertEqual(output, None)

  def test_run_fastml(self):
    builder = TreeBuilder()

    fasta_filename = os.path.join(test_data(), 'animals.mfa')
    tree_filename = os.path.join(test_data(), 'animals.terminal_nodes.newick')
    output_directory = tempfile.mkdtemp()

    # fastml doesn't like comments in fasta files, guesses it is a different
    # format and gets confused when a sequence has a '>' in it.  I'm therefore
    # creating a temporary file without the comments in it.
    fasta_file_without_comments = create_fasta_without_comments(fasta_filename)

    fastml_stdout, fastml_stderr = builder._run_fastml('fastml', {},
                                                       tree_filename,
                                                       fasta_file_without_comments,
                                                       output_directory)

    self.assertTrue(os.path.isfile(os.path.join(output_directory,
                                                'all_nodes.newick')))
    self.assertTrue(os.path.isfile(os.path.join(output_directory,
                                                'all_nodes.mfa')))

    shutil.rmtree(output_directory)
    os.remove(fasta_file_without_comments)

  @patch('biopericles.TreeBuilder.tempfile')
  def test_add_hereditary_nodes(self, temp_mock):
    builder = TreeBuilder()

    output_directory = tempfile.mkdtemp()
    temp_mock.mkdtemp.return_value = output_directory
    temp_mock.NamedTemporaryFile.side_effect = \
        (tempfile.NamedTemporaryFile(dir=output_directory, delete=False) for i in xrange(10))

    fasta_filename = os.path.join(test_data(), 'animals.mfa')
    with open(fasta_filename) as fasta_file:
      builder.load_fasta_sequences(fasta_file)

    tree_filename = os.path.join(test_data(), 'animals.terminal_nodes.newick')
    builder.tree = Bio.Phylo.read(tree_filename, 'newick')

    builder.add_hereditary_nodes()

    self.assertIsInstance(builder.tree, Bio.Phylo.Newick.Tree)
    tree_nodes = sorted([node.name for node in builder.tree.find_clades()])

    self.assertItemsEqual(tree_nodes, builder.sequences.keys())
    self.assertEqual(len(tree_nodes), 19)
    self.assertFalse(os.path.isdir(output_directory))

  def test_write_output(self):
    builder = TreeBuilder()
    builder.sequences = {"seq1": "a sequence object"}
    builder.tree = "Tree"
    builder.sequences_output_file = StringIO()
    builder.tree_output_file = StringIO()

    builder._write_tree = MagicMock()
    builder._write_sequences = MagicMock()

    builder.write_output()

    builder._write_tree.assert_called_with("Tree",
                                           builder.tree_output_file)
    builder._write_sequences.assert_called_with(["a sequence object"],
                                                builder.sequences_output_file)

  def test_write_tree(self):
    builder = TreeBuilder()

    tree_filename = os.path.join(test_data(), 'animals.terminal_nodes.newick')
    tree = Bio.Phylo.read(tree_filename, 'newick')

    output = StringIO()

    builder._write_tree(tree, output)

    with open(tree_filename, 'r') as tree_file:
      expected = tree_file.read()
    output.seek(0)

    self.assertEqual(output.read(), expected)

  def test_write_sequences(self):
    builder = TreeBuilder()

    fasta_filename = os.path.join(test_data(), 'animals.mfa')
    fasta_file = open(fasta_filename, 'r')
    sequences = Bio.SeqIO.parse(fasta_file, 'fasta')

    output = StringIO()

    builder._write_sequences(sequences, output)

    fasta_file.close()

    with open(fasta_filename, 'r') as fasta_file:
      expected = fasta_file.readlines()
      expected = "".join(expected[5:]) # Remove comments from the top
    output.seek(0)

    self.assertEqual(output.read(), expected)
