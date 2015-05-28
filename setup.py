import os
from setuptools import setup, find_packages
import multiprocessing

def readme():
  with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as f:
    return f.read()

setup(name='biopericles',
      version='0.1.0',
      description='Extracts features from genome sequences and classifies them into specified groupings',
      long_description=readme(),
      url='https://github.com/sanger-pathogens/BioPericles',
      author='Ben Taylor, Andrew J. Page',
      author_email='ben.taylor@sanger.ac.uk, ap13@sanger.ac.uk',
      scripts=['scripts/calculate-ancestral-sequences',
               'scripts/cluster-sequences',
               'scripts/get-cluster-consensus',
               'scripts/get-snp-features',
               'scripts/train-classifier']
      include_package_data=True,
      install_requires=[
        'biopython',
        'mock',
        'numpy',
        'PyVCF',
        'scikit-learn',
        'scipy'
      ],
      test_suite='nose.collector',
      tests_require=[
        'nose',
        'mock'
      ],
      license='GPLv3',
      packages=find_packages(),
      classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
      ]
)
