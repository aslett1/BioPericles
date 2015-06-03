#!/bin/bash

set -x
set -e

start_dir=$(pwd)

RAXML_VERSION="8.1.20"
FASTML_VERSION="3.1"
SNP_SITES_VERSION="1.5.0"
GENOMETOOLS_VERSION="1.5.5"

RAXML_DOWNLOAD_URL="https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz"
FASTML_DOWNLOAD_URL="http://fastml.tau.ac.il/source/FastML.v${FASTML_VERSION}.tgz"
SNP_SITES_DOWNLOAD_URL="https://github.com/sanger-pathogens/snp_sites/archive/${SNP_SITES_VERSION}.tar.gz"
GENOMETOOLS_DOWNLOAD_URL="https://github.com/genometools/genometools/archive/v${GENOMETOOLS_VERSION}.tar.gz"

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}

download $RAXML_DOWNLOAD_URL "raxml-${RAXML_VERSION}.tgz"
download $FASTML_DOWNLOAD_URL "fastml-${FASTML_VERSION}.tgz"
download $SNP_SITES_DOWNLOAD_URL "snp_sites-${SNP_SITES_VERSION}.tgz"
download $GENOMETOOLS_DOWNLOAD_URL "genometools-${GENOMETOOLS_VERSION}.tgz"

# Update dependencies
if [ "$TRAVIS" = 'true' ]; then
  echo "Using Travis's apt plugin"
else
  sudo apt-get update -q
  sudo apt-get install -y -q g++ \
  	                   zlib1g-dev check autoconf libtool git \
  		           python-dev libblas-dev liblapack-dev gfortran
fi

# Build all the things
cd $build_dir

## RAxML
raxml_dir=$(pwd)/"standard-RAxML-${RAXML_VERSION}"
if [ ! -d $raxml_dir ]; then
  tar xzf raxml-${RAXML_VERSION}.tgz
fi
cd $raxml_dir
if [ -e "${raxml_dir}/raxmlHPC" ]; then
  echo "Already build RAxML; skipping build"
else
  make -f Makefile.gcc
fi

cd $build_dir

## FASTML
fastml_dir=$(pwd)/"FastML.v${FASTML_VERSION}"
if [ ! -d $fastml_dir ]; then
  tar xzf fastml-${FASTML_VERSION}.tgz
fi
cd $fastml_dir
if [ -e "${fastml_dir}/programs/fastml/fastml" ]; then
  echo "Already build FASTML; skipping build"
else
  make
fi

cd $build_dir

## SNP SITES
snp_sites_dir=$(pwd)/snp_sites-${SNP_SITES_VERSION}
if [ ! -d $snp_sites_dir ]; then
  tar xzf snp_sites-${SNP_SITES_VERSION}.tgz
fi
cd $snp_sites_dir
if [ -e "${snp_sites_dir}/src/.libs/snp-sites" ]; then
  echo "Already build SNP Sites; skipping build"
else
  if [ ! -d m4 ]; then mkdir m4; fi
  autoreconf -i
  ./configure
  make
fi
if [ ! -e "${snp_sites_dir}/src/.libs/snp_sites" ]; then
  ln -s ${snp_sites_dir}/src/.libs/snp-sites ${snp_sites_dir}/src/.libs/snp_sites
fi

cd $build_dir

## GenomeTools
genometools_dir=$(pwd)/genometools-${GENOMETOOLS_VERSION}
if [ ! -d $genometools_dir ]; then
  tar xzf genometools-${GENOMETOOLS_VERSION}.tgz
fi
cd $genometools_dir
if [ -e "${genometools_dir}/bin/gt" ]; then
  echo "Already built GenomeTools; skipping build"
else
  make cairo=no
fi

cd $build_dir

# Setup environment variables
update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${raxml_dir}
update_path ${fastml_dir}/programs/fastml
export LD_LIBRARY_PATH=${snp_sites_dir}/src/.libs
update_path ${snp_sites_dir}/src/.libs
update_path ${genometools_dir}/bin

cd $start_dir

set +x
set +e
