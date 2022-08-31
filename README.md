# ipaPy2
Python implementation of the Integrated Probabilistic Annotation (IPA) - A Bayesian annotation method for LC/MS data integrating biochemical relations,
isotope patterns and adduct formation.

***
## Compiling from source (macOS)

1. create a folder in which you want to put the library
```
mkdir IPA
cd IPA
```
2. download the library
```
sudo apt install git
git clone https://github.com/francescodc87/ipaPy2
cd ipaPy2
```
3. create and activate a virtual environment for your folder and install the necessary libraries
```
python3 -m venv ipaPy2
source ipaPy2/bin/activate
pip install wheel
pip install setuptools
pip install twine
pip install pytest==4.4.1
pip install pytest-runner==4.4
```
4. run tests (optional)
```
python setup.py pytest
```
5. build your library
```
python setup.py bdist_wheel
```
6. The wheel file will be stored in the \dist folder. You can install the library in a new terminal as follows
```
pip install /path/to/wheelfile.whl
```

