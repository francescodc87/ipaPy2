# ipaPy2
Python implementation of the Integrated Probabilistic Annotation (IPA) - A Bayesian annotation method for LC/MS data integrating biochemical relations,
isotope patterns and adduct formation.

## Installation
ipaPy2 requires Python 3.9 or higher

### Install via pip (recommended )
NOT THERE YET!

### Compiling from source (macOS)
1. create a folder in which you want to put the library
```
mkdir IPA
cd IPA
```
2. download the library. If Homebrew is not installed in your machine, you can install it from here https://brew.sh 
```
brew install git
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

### Compiling from source (Linux)
1. create a folder in which you want to put the library
```
mkdir IPA
cd IPA
```
2. download the library
```
sudo apt-get install git
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

### Compiling from source (Windows)
to be added

## Database preparation
At the following link you can find a ready to use database for the ipaPy2 package. The fragmentaiton data has been collected from the MoNa database (https://mona.fiehnlab.ucdavis.edu), and it only includes the spectra acquired with a Q-exactive instrument.
One of the most powerful fetures of the IPA method is that it is able to integrate the knowledge gained from previous experiments in the anntoation process. Such knowledge must be included in the database and the following Jupyter Notebook detailes how the infomation can be added to the provided database template.
<br />
[Database preparation](tutorials/database_preparation.ipynb)

## Data preparation
In order to properly use this package, the processed data coming from an untargeted metabolomics experiment must be properly prepared.
The following tutorial details the proper formatting for the data.
<br />
[Dataset preparation](tutorials/dataset_preparation.ipynb)

## Usage
The Integrated Probabilistic Annotation (IPA) method can be applied in different situations. Below you can find a Jupyter notebook tutorial
showing how to apply the IPA method in different scenarios.
In order to run the tutorial, the following example dataset should be downloaded and the Jupyter notebooks should be run within the same folder.
https://link-to-dataset.com
<br />
[Usage tutorial](tutorials/Usage_tutorial.ipynb)


