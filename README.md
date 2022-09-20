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
```
```
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

## Databases
One of the most powerful fetures of the IPA method is that it is able to integrate the knowledge gained from previous experiments in the anntoation process. There are three files that are used as database:

**1. adducts file (required)**
<br />
The ipaPy2 library requires a file contains all the information required for the computation of the adducts. An adducts.csv file is provided with the package [here](DB/adducts.csv). The file contains the most common adducts. If any exotic adduct (or in-source fragment) needs to be considered, the user must modify the file accordingly. The format required for the adducts file is shown below. 



```python
import pandas as pd
import numpy as np
adducts = pd.read_csv('DB/adducts.csv')
adducts.head()
```



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
- **adductsPos**: list of adducts that should be considered in positive mode for this entry (e.g.,'M+Na;M+H;M+') - *necessary*
- **adductsNeg**: list of adducts that should be considered in Negative mode for this entry (e.g.,'M-H;M-2H') - *necessary*
- **description**: comments on the entry - *optional*
- **pk**: previous knowledge on the likelihood of this compoud to be present in the sample analyse. The value has to be between 1 (compound likely to be present in the sample) and 0 (compound cannot be present in the sample).
- **MS2**: id for the MS2 database entries related to this compound - *optional*
- **reactions**: list of reactions ids involving this compound (e.g., 'R00010 R00015 R00028'). If required, these can be used to find possible biochemical connections - *optional* 

The column names must be the ones reported here.
While the users are strongly advised to build their own *ad-hoc* database, [here](DB/IPA_MS1.csv) you can find a relatively big example database.


```python
DB = pd.read_csv('DB/IPA_MS1.csv')
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
      <th>adductsPos</th>
      <th>adductsNeg</th>
      <th>description</th>
      <th>pk</th>
      <th>MS2</th>
      <th>reactions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>C00002</td>
      <td>ATP</td>
      <td>C10H16N5O13P3</td>
      <td>InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
      <td>NaN</td>
      <td>1</td>
      <td>EMBL-MCF_spec365637_1</td>
      <td>R00002 R00076 R00085 R00086 R00087 R00088 R000...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>C00003</td>
      <td>NAD+</td>
      <td>C21H28N7O14P2</td>
      <td>InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
      <td>NaN</td>
      <td>1</td>
      <td>EMBL-MCF_specxxxxx_10</td>
      <td>R00023 R00090 R00091 R00092 R00093 R00094 R000...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C00004</td>
      <td>NADH</td>
      <td>C21H29N7O14P2</td>
      <td>InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
      <td>NaN</td>
      <td>1</td>
      <td>NaN</td>
      <td>R00023 R00090 R00091 R00092 R00093 R00094 R000...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>C00005</td>
      <td>NADPH</td>
      <td>C21H30N7O17P3</td>
      <td>InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
      <td>NaN</td>
      <td>1</td>
      <td>NaN</td>
      <td>R00105 R00106 R00107 R00108 R00109 R00111 R001...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C00006</td>
      <td>NADP+</td>
      <td>C21H29N7O17P3</td>
      <td>InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
      <td>NaN</td>
      <td>1</td>
      <td>EMBL-MCF_specxxxxxx_45</td>
      <td>R00104 R00106 R00107 R00108 R00109 R00111 R001...</td>
    </tr>
  </tbody>
</table>
</div>



This example databases was obtained considering the [KEGG database](https://www.genome.jp/kegg/compound/), the [Natural Products Atlas database](https://www.npatlas.org) and the [MoNa database](https://mona.fiehnlab.ucdavis.edu) (only compounds having at least one fragmentation spectra obtained with a Qexactive).
For each entry, only a handful of the most common adducts is considered.
To fully exploit the IPA method, it is strongly recommended to constantly update the database when new knowledge is gained from previous experience. Providing a retention time window for compounds previously detected with the analytical system at hand it is particularly useful.
For the sake of the example in this tutorial, a reduced example database is also provided.


```python
DB = pd.read_csv('DB/DB_test_pos.csv')
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
      <th>adductsPos</th>
      <th>adductsNeg</th>
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
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
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
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
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
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
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
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
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
      <td>M+H;M+Na;M+2H;2M+H</td>
      <td>M-H;2M-H;M-2H;3M-H</td>
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
This new implementation of the IPA method also allows the user to include MS2 data in the annotation pipeline.
In order to exploit this functionality a MS2 spectra database must be provided.
The MS2 database must be provided as a pandas dataframe including the following columns in this exact order:
- **compound_id**: unique id for each compound, it must match with the ids used in the MS1 database - *necessary*
- **id**: Unique id for the single entry (i.e., spectra) of the databse *necessary*
- **name**: compund name (e.g., 'D-Glucose') - *necessary*
- **formula**: chemical formula (e.g., 'C6H12O6') - *necessary*
- **inchi**: inchi string - *optional*
- **precursorType**: the adduct form of the precursor ion (e.g., 'M+H') - *necessary*
- **instrument**: the type of instrument the spectrum was acquired with - *optional*
- **collision.energy**: the collision energy level used to acquire the spectrum (e.g., '15') - *necessary*
- **spectrum**: The actual spectrum in the form of a string in the following format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'

It is necessary that the user uses a MS2 database specific to the instrument used to acquire the data.
The MS2 database found [here](https://drive.google.com/file/d/15qduvtE8aSAAUCf1FE4ojcVLaTw-B2W6/view?usp=sharing), contains all the MS2 spectra found in the [MoNa](https://mona.fiehnlab.ucdavis.edu) database acquired with a Qexactive. This is a relatively big file, and for the sake of this tutorial a drastically reduced version of it has been included within this repository, and can be foudn [here](DB/DBMS2_test_pos.csv).



```python
DBMS2 = pd.read_csv('DB/DBMS2_test_pos.csv')
DBMS2.head()
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
      <th>compound_id</th>
      <th>id</th>
      <th>name</th>
      <th>formula</th>
      <th>inchi</th>
      <th>precursorType</th>
      <th>instrument</th>
      <th>collision.energy</th>
      <th>spectrum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec103039</td>
      <td>L-valine</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>55.0550575256:5.821211 57.0581207275:0.385600 ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec353465</td>
      <td>L-valine</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028053042:0.000000 49.5031356971:0.000000 ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>EMBL-MCF_specxxxxxx_11</td>
      <td>EMBL-MCF_spec27828</td>
      <td>L-valine</td>
      <td>C5H11NO2</td>
      <td>InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1...</td>
      <td>M-H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>55.0173454285:0.819357 58.0282363892:0.155430 ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec96902</td>
      <td>L-proline</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>35</td>
      <td>50.5765228271:0.040013 51.3066940308:0.039949 ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>EMBL-MCF_specxxxxx_7</td>
      <td>EMBL-MCF_spec353568</td>
      <td>L-proline</td>
      <td>C5H9NO2</td>
      <td>InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...</td>
      <td>M+H</td>
      <td>Thermo Q-Exactive Plus</td>
      <td>30</td>
      <td>49.5028215674:0.000000 49.5031519602:0.000000 ...</td>
    </tr>
  </tbody>
</table>
</div>



## Data preparation
Before using the ipaPy2 package, the processed data coming from an untargeted metabolomics experiment must be properly prepared.

**1. MS1 data**
The data must be organized in a pandas dataframe containing the following columns:
- **ids**: an unique numeric id for each mass spectrometry feature feature
- **rel.ids**: relation ids. Features must be clustered based on correlation/peak shape/retention time. Features in the same cluster are likely to come from the same metabolite.
- **mzs**: mass-to-charge ratios, usually the average across different samples.
- **RTs**: retention times in seconds, usually the average across different samples.
- **Int**: epresentative (e.g., maximum or average) intensity detected for each feature across samples (either peak area or peak intensity)


Below is reported an example:


```python
df=pd.read_csv('example/df_test_pos.csv')
df.head()
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
      <th>ids</th>
      <th>rel.ids</th>
      <th>mzs</th>
      <th>RTs</th>
      <th>Int</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0</td>
      <td>116.070544</td>
      <td>45.770423</td>
      <td>2.170017e+09</td>
    </tr>
    <tr>
      <th>1</th>
      <td>88</td>
      <td>0</td>
      <td>117.073678</td>
      <td>45.787586</td>
      <td>1.256520e+08</td>
    </tr>
    <tr>
      <th>2</th>
      <td>501</td>
      <td>0</td>
      <td>231.133673</td>
      <td>46.183948</td>
      <td>2.519223e+07</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4429</td>
      <td>0</td>
      <td>232.136923</td>
      <td>46.176715</td>
      <td>2.635594e+06</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2</td>
      <td>1</td>
      <td>104.106830</td>
      <td>40.843309</td>
      <td>1.889172e+09</td>
    </tr>
  </tbody>
</table>
</div>



The clustering of the features is a necessary and must be performed before running the IPA method. For this step, the use of widely used data processing software such as [mzMatch](https://github.com/UoMMIB/mzmatch.R) and [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html) is recommended.
Nevertheless, the ipaPy2 library provides a function able to perform such step, starting from a dataframe containing the measured intensities across several samples (at least 3 samples, the more samples the better). 
Such dataframe should be organized as follows:


```python
df2=pd.read_csv('example/df_test_pos_not_clustered.csv')
df2.head()
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
      <th>ids</th>
      <th>mzs</th>
      <th>RTs</th>
      <th>sample1</th>
      <th>sample2</th>
      <th>sample3</th>
      <th>sample4</th>
      <th>sample5</th>
      <th>sample6</th>
      <th>sample7</th>
      <th>sample8</th>
      <th>sample9</th>
      <th>sample10</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>116.070544</td>
      <td>45.770423</td>
      <td>1.003660e+09</td>
      <td>1.299828e+09</td>
      <td>1.878029e+09</td>
      <td>1.778238e+09</td>
      <td>1.715394e+09</td>
      <td>4.340034e+08</td>
      <td>1.586635e+09</td>
      <td>2.170017e+09</td>
      <td>1.312151e+09</td>
      <td>2.051875e+09</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>104.106830</td>
      <td>40.843309</td>
      <td>3.778343e+08</td>
      <td>8.721901e+08</td>
      <td>8.353805e+08</td>
      <td>1.889172e+09</td>
      <td>1.114844e+09</td>
      <td>1.296362e+09</td>
      <td>7.361379e+08</td>
      <td>7.386887e+08</td>
      <td>9.546864e+08</td>
      <td>6.969054e+08</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>118.085998</td>
      <td>43.584638</td>
      <td>5.984715e+08</td>
      <td>1.399106e+09</td>
      <td>2.831220e+08</td>
      <td>1.415610e+09</td>
      <td>7.557607e+08</td>
      <td>7.800359e+08</td>
      <td>8.949854e+08</td>
      <td>5.074069e+08</td>
      <td>6.854525e+08</td>
      <td>1.000501e+09</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>166.086047</td>
      <td>143.321396</td>
      <td>1.390905e+09</td>
      <td>1.047887e+09</td>
      <td>1.053413e+09</td>
      <td>2.781809e+08</td>
      <td>1.037486e+09</td>
      <td>1.117700e+09</td>
      <td>6.153332e+08</td>
      <td>1.215932e+09</td>
      <td>1.264092e+09</td>
      <td>1.370995e+09</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>132.101745</td>
      <td>89.387202</td>
      <td>6.071912e+08</td>
      <td>1.014152e+09</td>
      <td>1.270735e+09</td>
      <td>1.069765e+09</td>
      <td>4.925938e+08</td>
      <td>4.087633e+08</td>
      <td>3.777945e+08</td>
      <td>2.541470e+08</td>
      <td>8.025257e+08</td>
      <td>3.544281e+08</td>
    </tr>
  </tbody>
</table>
</div>



clusting function description


```python
#clustering function usage
```

**2. MS2 data**
If fragmentation data was acquired during the experiment, it can be included in the IPA annotation process.
To do so, the data must be organized in a pandas dataframe containing the following columns:
- **id**: an unique id for each feature for which the MS2 spectrum was acquired (same as in MS1)
- **spectrum**: string containing the spectrum inforamtion in the following format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'
- **ev**: collision energy used to aquire the fragmentation spectrum

Below is reported an example:


```python
dfMS2=pd.read_csv('example/MS2data_example.csv')
dfMS2.head()
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
      <th>spectrum</th>
      <th>ev</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>51.3066132836457:0.884272376680125 59.96532241...</td>
      <td>35</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>51.3066132836457:0.884272376680125 59.96532241...</td>
      <td>15</td>
    </tr>
    <tr>
      <th>2</th>
      <td>90</td>
      <td>62.4153253406374:0.743812036877455 63.93291389...</td>
      <td>35</td>
    </tr>
    <tr>
      <th>3</th>
      <td>992</td>
      <td>50.983321052233:0.973529955385613 53.039006800...</td>
      <td>35</td>
    </tr>
    <tr>
      <th>4</th>
      <td>3</td>
      <td>55.0551847656264:5.67780579195993 57.058126021...</td>
      <td>35</td>
    </tr>
  </tbody>
</table>
</div>



## Usage
The Integrated Probabilistic Annotation (IPA) method can be applied in different situations. Below you can find a Jupyter notebook tutorial
showing how to apply the IPA method in different scenarios.
In order to run the tutorial, the following example dataset should be downloaded and the Jupyter notebooks should be run within the same folder.
https://link-to-dataset.com
<br />
[Usage tutorial](tutorials/Usage_tutorial.ipynb)


