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

## Database
One of the most powerful fetures of the IPA method is that it is able to integrate the knowledge gained from previous experiments in the anntoation process. There are three files that are used as database:

**1. adducts file (required)**
<br />
This file contains all the information required for the computation of the adducts. An adducts.csv file is provided with the package [here](DB/adducts.csv). The file contains the most common adducts. If any exotic adduct (or in-source fragment) needs to be considered, the user must modify the file accordingly. The format required for the adducts file is shown below. 



```python
import pandas as pd
import numpy as np
adducts = pd.read_csv('DB/adducts.csv')
adducts.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Name</th>
      <th>calc</th>
      <th>Charge</th>
      <th>Mult</th>
      <th>Mass</th>
      <th>Ion_mode</th>
      <th>Formula_add</th>
      <th>Formula_ded</th>
      <th>Multi</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>M+H</td>
      <td>M+1.007276</td>
      <td>1</td>
      <td>1</td>
      <td>1.007276</td>
      <td>positive</td>
      <td>H1</td>
      <td>FALSE</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>M+NH4</td>
      <td>M+18.033823</td>
      <td>1</td>
      <td>1</td>
      <td>18.033823</td>
      <td>positive</td>
      <td>N1H4</td>
      <td>FALSE</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>M+Na</td>
      <td>M+22.989218</td>
      <td>1</td>
      <td>1</td>
      <td>22.989218</td>
      <td>positive</td>
      <td>Na1</td>
      <td>FALSE</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>M+K</td>
      <td>M+38.963158</td>
      <td>1</td>
      <td>1</td>
      <td>38.963158</td>
      <td>positive</td>
      <td>K1</td>
      <td>FALSE</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>M+</td>
      <td>M-0.00054858</td>
      <td>1</td>
      <td>1</td>
      <td>-0.000549</td>
      <td>positive</td>
      <td>FALSE</td>
      <td>FALSE</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



**2. MS1 database file (required)**
<br />
The IPA method requires a pandas dataframe containing the database against which the annotation is performed.
Such dataframe must contain the following columns in this exact order (optional colums can have empty fields):
- **id**: unique id of the database entry (e.g., 'C00031') - *necessary*
- **name**: compund name (e.g., 'D-Glucose') - *necessary*
- **formula**: chemical formula (e.g., 'C6H12O6') - *necessary*
- **inchi**: inchi string - *optional*
- **smiles**: smiles string - *optional*
- **RT**: if known, retention time range (in seconds) where this compound is expected to elute (e.g., '30;60') - *optional*
- **adducts**: list of adducts that should be considered for this entry (e.g.,'M+Na;M+H;M+')
- **description**: comments on the entry - *optional*
- **pk**: previous knowledge on the likelihood of this compoud to be present in the sample analyse. The value has to be between 1 (compound likely to be present in the sample) and 0 (compound cannot be present in the sample).
- **MS2**: id for the MS2 database entries related to this compound - *optional*
- **reactions**: list of reactions ids involving this compound (e.g., 'R00010 R00015 R00028')- *optional* 

The database is very important and it is strongly advised to update it in order to take full advantage of the IPA capabilities.
Here you can find an example DB...

A very small version of the DB with some information already filled in has been built only for the sake of the tutorial and can be found here...



```python
DB = pd.read_csv('DB/DB_test_pos.csv')
DB.replace(np.nan,None)
DB.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>name</th>
      <th>formula</th>
      <th>inchi</th>
      <th>smiles</th>
      <th>RT</th>
      <th>adducts</th>
      <th>description</th>
      <th>pk</th>
      <th>MS2</th>
      <th>reactions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>C00079</td>
      <td>L-Phenylalanine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-...</td>
      <td>NaN</td>
      <td>120;160</td>
      <td>M+H;M+2H;M+Na;2M+H;M+</td>
      <td>NaN</td>
      <td>1</td>
      <td>UA005501_1</td>
      <td>R00686 R00688 R00689 R00690 R00691 R00692 R006...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>C00082</td>
      <td>L-Tyrosine</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>NaN</td>
      <td>50;90</td>
      <td>M+H;M+2H;M+Na;2M+H;M+</td>
      <td>NaN</td>
      <td>1</td>
      <td>UA005601_1</td>
      <td>R00031 R00728 R00729 R00730 R00731 R00732 R007...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C00114</td>
      <td>Choline</td>
      <td>C5H14NO</td>
      <td>InChI=1S/C5H14NO/c1-6(2,3)4-5-7/h7H,4-5H2,1-3H...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+2H;M+Na;2M+H;M+</td>
      <td>NaN</td>
      <td>1</td>
      <td>NaN</td>
      <td>R01021 R01022 R01023 R01025 R01026 R01027 R010...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>C00123</td>
      <td>L-Leucine</td>
      <td>C6H13NO2</td>
      <td>InChI=1S/C6H13NO2/c1-4(2)3-5(7)6(8)9/h4-5H,3,7...</td>
      <td>NaN</td>
      <td>70;110</td>
      <td>M+H;M+2H;M+Na;2M+H;M+</td>
      <td>NaN</td>
      <td>1</td>
      <td>NaN</td>
      <td>R01088 R01089 R01090 R01091 R02552 R03657 R084...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C00148</td>
      <td>L-Proline</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>NaN</td>
      <td>35;55</td>
      <td>M+H;M+2H;M+Na;2M+H;M+</td>
      <td>NaN</td>
      <td>1</td>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>R00135 R00671 R01246 R01248 R01249 R01251 R012...</td>
    </tr>
  </tbody>
</table>
</div>



**3. MS2 database file (only required is MS2 data is available)**
<br />

At the following link you can find a ready to use database for the ipaPy2 package. The fragmentaiton data has been collected from the MoNa database (https://mona.fiehnlab.ucdavis.edu), and it only includes the spectra acquired with a Q-exactive instrument.
One of the most powerful fetures of the IPA method is that it is able to integrate the knowledge gained from previous experiments in the anntoation process. Such knowledge must be included in the database and the following Jupyter Notebook detailes how the infomation can be added to the provided database template.
<br />
[Database preparation](tutorials/database_preparation.ipynb)


```python
DBMS2 = pd.read_csv('DBMS2_test_pos.csv')
DBMS2
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>MonaID</th>
      <th>id</th>
      <th>Name</th>
      <th>formula</th>
      <th>inchi</th>
      <th>precursorType</th>
      <th>instrument</th>
      <th>collision.energy</th>
      <th>spectrum</th>
      <th>KEGG</th>
      <th>CAS</th>
      <th>PubChem</th>
      <th>ChEBI</th>
      <th>KNApSAcK</th>
      <th>LIPIDMAPS</th>
      <th>LipidBank</th>
      <th>chemspider</th>
      <th>HMDB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SM827301_1</td>
      <td>SM827301</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Plus Orbitrap Thermo Scientific</td>
      <td>35</td>
      <td>53.0388:0.282352 65.0388:0.272609 77.0386:0.90...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004765</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.810675:0.011222 51.102235:0.013754 51.50411...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004766</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.291872:0.286642 50.314117:0.201290 50.39813...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004767</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.077123:11.739013 50.739286:7.700284 50.8759...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004768</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+K</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.348437:0.012559 50.565674:0.012258 50.65556...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>5</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004769</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+K</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.175086:0.183738 50.179197:0.138681 50.33830...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>6</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004770</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+K</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.276741:0.255358 50.459141:0.356575 50.68504...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004771</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+NH4</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.022209:0.005153 50.390025:0.010307 51.00674...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004772</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+NH4</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.219187:0.035615 50.403890:0.053475 50.40557...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>9</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004773</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+NH4</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.015966:0.082126 50.123059:0.031558 50.13416...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004774</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.360543:0.016865 50.738460:0.034515 50.89488...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>11</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004775</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.501417:0.049480 50.650735:0.050883 50.81192...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>12</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004776</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.107105:0.064235 50.639921:0.043468 50.88637...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>13</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013240</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.209323:0.043050 50.548028:0.023113 50.58551...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>14</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013241</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.251420:0.114624 50.466517:0.053005 50.69270...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>15</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013242</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+H</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.000829:0.109223 50.601142:0.130240 51.01419...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>16</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013243</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.272464:0.024599 50.554612:0.012395 51.17923...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>17</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013244</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.621130:0.020986 50.866899:0.023192 51.31972...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013245</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M+Na</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.027019:0.027762 50.307668:0.016216 50.41594...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>19</th>
      <td>UA005501_1</td>
      <td>UA005501</td>
      <td>L-Phenylalanine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-...</td>
      <td>M+H</td>
      <td>Q Exactive Plus Orbitrap Thermo Scientific</td>
      <td>50</td>
      <td>79.0542:1.317033 91.0542:1.075051 93.0698:2.74...</td>
      <td>C00079</td>
      <td>63-91-2</td>
      <td>NaN</td>
      <td>58095.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>20</th>
      <td>UA005601_1</td>
      <td>UA005601</td>
      <td>L-Tyrosine</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M+H</td>
      <td>Q Exactive Plus Orbitrap Thermo Scientific</td>
      <td>50</td>
      <td>91.0542:49.584656 95.0491:19.730330 103.0542:1...</td>
      <td>D00022</td>
      <td>60-18-4</td>
      <td>NaN</td>
      <td>58315.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>21</th>
      <td>SM827301_1</td>
      <td>LU102301</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>15</td>
      <td>81.0699:0.674506 91.0542:0.318345 94.0413:0.31...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>22</th>
      <td>SM827301_1</td>
      <td>LU102302</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>30</td>
      <td>58.0288:0.129172 81.0698:0.749676 91.0542:0.23...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>23</th>
      <td>SM827301_1</td>
      <td>LU102303</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>45</td>
      <td>58.0288:0.550849 65.0386:0.248859 79.0543:0.58...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>24</th>
      <td>SM827301_1</td>
      <td>LU102304</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>60</td>
      <td>53.0386:0.191444 55.0178:0.192629 58.0288:1.40...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>25</th>
      <td>SM827301_1</td>
      <td>LU102305</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>75</td>
      <td>53.0023:0.143043 53.0386:0.895701 55.0178:0.88...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>26</th>
      <td>SM827301_1</td>
      <td>LU102306</td>
      <td>Benzocaine</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c1-2-12-9(11)7-3-5-8(10)6-4-...</td>
      <td>M+H</td>
      <td>Q Exactive Orbitrap (Thermo Scientific)</td>
      <td>90</td>
      <td>53.0021:0.652513 53.0386:2.650298 55.0178:1.79...</td>
      <td>D00552</td>
      <td>94-09-7</td>
      <td>NaN</td>
      <td>116735.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>27</th>
      <td>EMBL-MCF_spec3xxxxx_2</td>
      <td>EMBL-MCF_spec352816</td>
      <td>Isoleucine</td>
      <td>C6H13NO2</td>
      <td>InChI=1S/C6H13NO2/c1-3-4(2)5(7)6(8)9/h4-5H,3,7...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5027348542:0.000000 49.5030652471:0.000000 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>28</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec96902</td>
      <td>L-PROLINE</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>50.5765228271:0.040013 51.3066940308:0.039949 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>29</th>
      <td>UA005501_1</td>
      <td>EMBL-MCF_spec98214</td>
      <td>L-PHENYLALANINE</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>51.0238380432:0.039736 72.9652786255:0.054871 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>30</th>
      <td>UA005601_1</td>
      <td>EMBL-MCF_spec103377</td>
      <td>L-TYROSINE</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>50.0159416199:0.062644 56.3880615234:0.047060 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>31</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec103039</td>
      <td>L-VALINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>55.0550575256:5.821211 57.0581207275:0.385600 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>32</th>
      <td>EMBL-MCF_spec107xxx_1</td>
      <td>EMBL-MCF_spec107285</td>
      <td>BETAINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>52.0572242737:0.028822 57.2626113892:0.028542 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>33</th>
      <td>EMBL-MCF_spec107xxx_1</td>
      <td>EMBL-MCF_spec107314</td>
      <td>BETAINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3</td>
      <td>M+Na</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>52.6068687439:0.061225 53.0006103516:0.079595 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>34</th>
      <td>EMBL-MCF_specxxxxxx_35</td>
      <td>EMBL-MCF_spec111716</td>
      <td>L-NORVALINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-2-3-4(6)5(7)8/h4H,2-3,6H2...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>55.0550384521:0.200344 58.0659370422:0.711922 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>35</th>
      <td>EMBL-MCF_spec128577_1</td>
      <td>EMBL-MCF_spec128577</td>
      <td>5-AMINOPENTANOATE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c6-4-2-1-3-5(7)8/h1-4,6H2,(H...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>53.0394134521:0.336045 55.0550384521:29.021750...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>36</th>
      <td>UA005601_1</td>
      <td>EMBL-MCF_spec353490</td>
      <td>L-TYROSINE</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028380308:0.000000 49.5031684236:0.000000 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>37</th>
      <td>UA005501_1</td>
      <td>EMBL-MCF_spec352576</td>
      <td>L-PHENYLALANINE</td>
      <td>C9H11NO2</td>
      <td>InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028156518:0.000000 49.5031460446:0.000000 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>38</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec353465</td>
      <td>L-VALINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028053042:0.000000 49.5031356971:0.000000 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>39</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec353568</td>
      <td>L-PROLINE</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028215674:0.000000 49.5031519602:0.000000 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>40</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004762</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>51.074606:0.009458 51.895476:0.010647 52.15682...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>41</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004763</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.122187:0.036305 50.278440:0.029221 50.45575...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>42</th>
      <td>VF-NPL-QEHF00xxxx_74</td>
      <td>VF-NPL-QEHF004764</td>
      <td>geniposidic acid</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c17-3-6-1-2-7-8(14(22)23)5-...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.382877:0.075648 50.664436:0.059475 51.01725...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>43</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013237</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>35</td>
      <td>50.292765:0.017220 50.343621:0.007022 50.39343...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>44</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013238</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>45</td>
      <td>50.220643:0.011954 50.401338:0.020494 50.45846...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>45</th>
      <td>VF-NPL-QEHF0132xx_2</td>
      <td>VF-NPL-QEHF013239</td>
      <td>swertiamarin</td>
      <td>C16H22O10</td>
      <td>InChI=1S/C16H22O10/c1-2-7-14(24-6-8-13(21)23-4...</td>
      <td>M-H</td>
      <td>Thermo Q Exactive HF</td>
      <td>65</td>
      <td>50.274791:0.018123 50.393134:0.031226 50.85927...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>46</th>
      <td>EMBL-MCF_spec3xxxxx_2</td>
      <td>EMBL-MCF_spec370588</td>
      <td>Isoleucine</td>
      <td>C6H13NO2</td>
      <td>InChI=1S/C6H13NO2/c1-3-4(2)5(7)6(8)9/h4-5H,3,7...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>50.9206542969:30.130281 58.476108551:32.731387...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>47</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec24658</td>
      <td>L-PROLINE</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>59.0125083923:0.125338 66.0334091187:1.831800 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>48</th>
      <td>UA005601_1</td>
      <td>EMBL-MCF_spec26382</td>
      <td>L-TYROSINE</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M-H2O+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>53.3883934021:0.053113 63.8395233154:0.054732 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>49</th>
      <td>UA005601_1</td>
      <td>EMBL-MCF_spec28025</td>
      <td>L-TYROSINE</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>72.0076065063:22.592507 73.0110626221:0.236782...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>50</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec27828</td>
      <td>L-VALINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>55.0173454285:0.819357 58.0282363892:0.155430 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>51</th>
      <td>EMBL-MCF_spec36931_1</td>
      <td>EMBL-MCF_spec36931</td>
      <td>5-OXO-D-PROLINE</td>
      <td>C5H7NO3</td>
      <td>InChI=1S/C5H7NO3/c7-4-2-1-3(6-4)5(8)9/h3H,1-2H...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>50.1261749268:0.060748 59.012348175:0.090858 7...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>52</th>
      <td>EMBL-MCF_spec38257_1</td>
      <td>EMBL-MCF_spec38257</td>
      <td>5-OXO-L-PROLINE</td>
      <td>C5H7NO3</td>
      <td>InChI=1S/C5H7NO3/c7-4-2-1-3(6-4)5(8)9/h3H,1-2H...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>51.0553016663:0.142479 65.5978164673:0.120717 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>53</th>
      <td>EMBL-MCF_specxxxxxx_35</td>
      <td>EMBL-MCF_spec41318</td>
      <td>L-NORVALINE</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-2-3-4(6)5(7)8/h4H,2-3,6H2...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>60.5583457947:0.049289 60.9735183716:0.199877 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>54</th>
      <td>EMBL-MCF_spec59716_1</td>
      <td>EMBL-MCF_spec59716</td>
      <td>NORLEUCINE</td>
      <td>C6H13NO2</td>
      <td>InChI=1S/C6H13NO2/c1-2-3-4-5(7)6(8)9/h5H,2-4,7...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>60.9735603333:0.383644 71.0126419067:0.234627 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>55</th>
      <td>UA005601_1</td>
      <td>EMBL-MCF_spec391818</td>
      <td>L-TYROSINE</td>
      <td>C9H11NO3</td>
      <td>InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>50.539855957:0.574067 53.675239563:0.649962 54...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>56</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec388914</td>
      <td>L-PROLINE</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>66.0333938599:6.892720 68.8465042114:7.550322 ...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>



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


