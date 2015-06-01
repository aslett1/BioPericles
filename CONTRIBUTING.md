Please update the version and run the [tests](run_tests.sh) before contributing.

If you have root in a clean Ubuntu machine (e.g. a Docker container) you
can install all of the dependencies using:

```
cd BioPericles/prototype
sudo ./install_dependencies.sh
cd ..
pip install -r requirements.txt
python setup.py develop
```
