from invoke import task, Collection
from invoke import run as invoke_run
import contextlib
import os
import logging
import requests

logging.basicConfig(format='%(message)s', level=logging.INFO)

RAXML_VERSION = "8.1.20"
FASTML_VERSION = "3.1"
SNP_SITES_VERSION = "1.5.0"

RAXML_DOWNLOAD_URL = \
    "https://github.com/stamatak/standard-RAxML/archive/v{version}.tar.gz".format(version=RAXML_VERSION)
FASTML_DOWNLOAD_URL = \
    "http://fastml.tau.ac.il/source/FastML.v{version}.tgz".format(version=FASTML_VERSION)
SNP_SITES_DOWNLOAD_URL = \
    "https://github.com/sanger-pathogens/snp_sites/archive/{version}.tar.gz".format(version=SNP_SITES_VERSION)

build_dir = os.path.join(os.getcwd(), "invoke_install_directory")

ns = Collection()

class RunException(Exception):
  def __init__(self, message, result):
    super(Exception, self).__init__(message)
    self.result = result

class DownloadException(Exception):
  def __init__(self, message, response):
    super(Exception, self).__init__(message)
    self.response = response

@contextlib.contextmanager
def cd(path):
  old_dir = os.path.abspath(os.getcwd())
  try:
    if os.path.abspath(path) != old_dir:
      logging.info("--> cd into %s" % path)
    os.chdir(path)
    yield
  finally:
    os.chdir(old_dir)

def run(command, *args, **kwargs):
  logging.info("--> {directory}$ {command}".format(directory=os.getcwd(),
                                                command=command))
  result = invoke_run(command, *args, **kwargs)
  logging.info(result.stdout)
  if result.return_code != 0:
    message = "Issue runnung {command} in {directory}.  Exited with {status}".format(command=command, 
                                                                                     directory=os.getcwd(),
                                                                                     status=result.return_code)
    raise RunException(message, result)

def download(url, path):
  logging.info("--> downloading {url} to {path}".format(url=url,
                                                        path=os.path.abspath(path)))
  with open(path, 'w') as output_file:
    response = requests.get(url, stream=True)
    if not response.ok:
      raise DownloadException("Could not download %s. Status code: %s" % (url, response.status_code),
                              response)
    for block in response.iter_content(1024):
      if not block:
        break
      output_file.write(block)

def say(message):
  logging.info("--> %s" % message)

def isinstalled(package):
  try:
    result = invoke_run("dpkg -s %s 2>/dev/null" % package)
    return result.return_code == 0
  except:
    return False

def apt_install(packages):
  if isinstance(packages, str):
    packages = packages.split()
  for package in packages:
    if not isinstalled(package):
      run("sudo apt-get install -y %s" % package)

def skip_task(task):
  logging.info("--> Task '%s' already complete, skipping" % task)

@task
def create_install_folder():
  if os.path.isdir(build_dir):
    skip_task("create install folder")
  else:
    run("mkdir %s" % build_dir)

@task
def delete_install_folder():
  if not os.path.isdir(build_dir):
    skip_task("delete install folder")
  else:
    run("rm -r %s" % build_dir)

@task(pre=[create_install_folder])
def download_raxml():
  with cd(build_dir):
    output_file = "raxml-{version}.tgz".format(version=RAXML_VERSION)
    if os.path.isfile(output_file):
      skip_task("download %s" % output_file)
    else:
      download(RAXML_DOWNLOAD_URL, output_file)

@task(pre=[create_install_folder])
def download_fastml():
  with cd(build_dir):
    output_file = "fastml-{version}.tgz".format(version=FASTML_VERSION)
    if os.path.isfile(output_file):
      skip_task("download %s" % output_file)
    else:
      download(FASTML_DOWNLOAD_URL, output_file)

@task(pre=[create_install_folder])
def download_snp_sites():
  with cd(build_dir):
    output_file = "snp_sites-{version}.tgz".format(version=SNP_SITES_VERSION)
    if os.path.isfile(output_file):
      skip_task("download %s" % output_file)
    else:
      download(SNP_SITES_DOWNLOAD_URL, output_file)

@task(pre=[download_raxml, download_fastml, download_snp_sites])
def download_source():
  pass

download_tasks = Collection('download')
download_tasks.add_task(download_source, 'all', default=True)
download_tasks.add_task(download_raxml, 'raxml')
download_tasks.add_task(download_fastml, 'fastml')
download_tasks.add_task(download_snp_sites, 'snp_sites')
ns.add_collection(download_tasks)

@task
def apt_update():
  run("sudo apt-get update")

@task(pre=[download_raxml])
def extract_raxml():
  with cd(build_dir):
    if os.path.isdir("standard-RAxML-%s" % RAXML_VERSION):
      skip_task("extract raxml")
    else:
      run("tar xzf raxml-%s.tgz" % RAXML_VERSION)

@task(pre=[extract_raxml])
def build_raxml():
  raxml_directory=os.path.join(build_dir, "standard-RAxML-%s" % RAXML_VERSION)
  with cd(raxml_directory):
    if os.path.isfile('raxmlHPC'):
      skip_task('build raxml')
    else:
      run("make -f Makefile.gcc")

@task
def install_raxml_dependencies():
  pass

@task(pre=[install_raxml_dependencies, build_raxml])
def install_raxml():
  raxml_directory=os.path.join(build_dir, "standard-RAxML-%s" % RAXML_VERSION)
  with cd(raxml_directory):
    run("sudo cp raxmlHPC /usr/local/bin/")

@task(pre=[apt_update])
def install_fastml_dependencies():
  apt_install("g++")

@task(pre=[download_fastml])
def extract_fastml():
  with cd(build_dir):
    if os.path.isdir("FastML.v%s" % FASTML_VERSION):
      skip_task("extract fastml")
    else:
      run("tar xzf fastml-%s.tgz" % FASTML_VERSION)

@task(pre=[extract_fastml])
def build_fastml():
  fastml_directory=os.path.join(build_dir, "FastML.v%s" % FASTML_VERSION)
  with cd(fastml_directory):
    if os.path.isfile(os.path.join('programs', 'fastml', 'fastml')):
      skip_task("build fastml")
    else:
      run('make')

@task(pre=[install_fastml_dependencies, build_fastml])
def install_fastml():
  fastml_directory=os.path.join(build_dir, "FastML.v%s" % FASTML_VERSION)
  with cd(fastml_directory):
    run("sudo cp %s /usr/local/bin" % os.path.join('programs', 'fastml', 'fastml'))

@task(pre=[apt_update])
def install_snp_sites_dependencies():
  apt_install('zlib1g-dev check autoconf libtool git')

@task(pre=[install_snp_sites_dependencies, download_snp_sites])
def build_snp_sites():
  snp_sites_directory = os.path.join(build_dir, "snp_sites")
  with cd(snp_sites_directory):
    if os.path.isfile(os.path.join("src", ".libs", "snp-sites")):
      skip_task("build snp_sites")
    else:
      if not os.path.isdir("m4"):
        run("mkdir m4")
      run("autoreconf -i")
      run("./configure")
      run("make")

@task(pre=[build_snp_sites])
def install_snp_sites():
  snp_sites_directory = os.path.join(build_dir, "snp_sites")
  with cd(snp_sites_directory):
    run("sudo make install")
    say("please run `export LD_LIBRARY_PATH=/usr/local/lib` in your terminal")

@task(pre=[install_raxml, install_fastml, install_snp_sites],
      post=[delete_install_folder])
def install_all():
  pass

install_tasks = Collection('install')
install_tasks.add_task(install_all, 'all', default=True)
install_tasks.add_task(install_raxml, 'raxml')
install_tasks.add_task(install_fastml, 'fastml')
install_tasks.add_task(install_snp_sites, 'snp_sites')
ns.add_collection(install_tasks)
