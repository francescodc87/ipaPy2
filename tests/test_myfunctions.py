from ipaPy2 import MS2compare
from ipaPy2 import ipa
import pandas as pd
import numpy as np
import pickle


def test_clusterFeatures():
    file = open('test_clusterFeatures.pkl', 'rb')
    df= pickle.load(file)
    expected = pickle.load(file)
    file.close()
    df2 = ipa.clusterFeatures(df)
    assert(df2.equals(expected))\


def test_map_isotope_patterns():
    file = open('test_map_isotope_patterns.pkl', 'rb')
    df= pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    df_out =df.copy()
    ipa.map_isotope_patterns(df_out,ionisation=1)
    assert(df_out.equals(expected1))
    
    df_out =df.copy()
    ipa.map_isotope_patterns(df_out,ionisation=-1)
    assert(df_out.equals(expected2))



def test_compute_all_adducts():
    file = open('test_compute_all_adducts.pkl', 'rb')
    adducts=pickle.load(file)
    DB=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    
    out = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=1)
    assert(out.equals(expected1))

    out = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=2)
    assert(out.equals(expected1))
    
    out = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=1)
    assert(out.equals(expected2))

    out = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=2)
    assert(out.equals(expected2))


def test_MS1annotation():
    file = open('test_MS1annotation.pkl', 'rb')
    df=pickle.load(file)
    allAdds=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    annotations=ipa.MS1annotation(df,allAdds,ppm=5,
                                  me = 5.48579909065e-04,ratiosd=0.9,
                                  ppmunk=None,ratiounk=None,ppmthr=None,
                                  pRTNone=None,pRTout=None,ncores=1)
    
    assert(annotations[1].equals(expected1[1]))
    assert(annotations[501].equals(expected1[501]))
    assert(annotations[4].equals(expected1[4]))
    assert(annotations[999].equals(expected1[999]))

    annotations=ipa.MS1annotation(df,allAdds,ppm=5,
                                  me = 5.48579909065e-04,ratiosd=0.9,
                                  ppmunk=None,ratiounk=None,ppmthr=None,
                                  pRTNone=None,pRTout=None,ncores=2)
    assert(annotations[1].equals(expected2[1]))
    assert(annotations[501].equals(expected2[501]))
    assert(annotations[4].equals(expected2[4]))
    assert(annotations[999].equals(expected2[999]))
    

def test_MSMSannotation():
    file = open('test_MS2annotation.pkl', 'rb')
    df=pickle.load(file)
    allAdds=pickle.load(file)
    dfMS2=pickle.load(file)
    DBMS2=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    annotations=ipa.MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm=5,
                               me=0.000548579909065,ratiosd=0.9,
                               ppmunk=None,ratiounk=None,ppmthr=None,
                               pRTNone=None,pRTout=None,mzdCS=0,ppmCS=10,
                               CSunk=0.7,evfilt=False,ncores=1)
    
    assert(annotations[1].equals(expected1[1]))
    assert(annotations[501].equals(expected1[501]))
    assert(annotations[4].equals(expected1[4]))
    assert(annotations[999].equals(expected1[999]))

    annotations=ipa.MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm=5,
                               me=0.000548579909065,ratiosd=0.9,
                               ppmunk=None,ratiounk=None,ppmthr=None,
                               pRTNone=None,pRTout=None,mzdCS=0,ppmCS=10,
                               CSunk=0.7,evfilt=False,ncores=2)
    
    assert(annotations[1].equals(expected1[1]))
    assert(annotations[501].equals(expected1[501]))
    assert(annotations[4].equals(expected1[4]))
    assert(annotations[999].equals(expected1[999]))


def test_Compute_Bio():
    file = open('test_Compute_Bio.pkl', 'rb')
    DB=pickle.load(file)
    annotations=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    Bio=ipa.Compute_Bio(DB,annotations,mode='reactions',ncores=1)
    assert(Bio.equals(expected1))
    Bio=ipa.Compute_Bio(DB,annotations,mode='reactions',ncores=2)
    assert(Bio.equals(expected1))
    Bio=ipa.Compute_Bio(DB,annotations,mode='connections',connections=['C4H2'],ncores=1)
    assert(Bio.equals(expected2))
    Bio=ipa.Compute_Bio(DB,annotations,mode='connections',connections=['C4H2'],ncores=2)
    assert(Bio.equals(expected2))
    
    
def test_Gibbs_sampler_add():
    file = open('test_Gibbs_sampler_add.pkl', 'rb')
    df=pickle.load(file)
    annotations=pickle.load(file)
    expected1 = pickle.load(file)
    file.close()
    
    ipa.Gibbs_sampler_add(df, annotations, noits=1000, burn=None, delta_add=.5, all_out=False, zs=None)
    
    ex1=expected1[1]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[501]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[4]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[999]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))

def test_Gibbs_sampler_bio():
    file = open('test_Gibbs_sampler_bio.pkl', 'rb')
    df=pickle.load(file)
    Bio=pickle.load(file)
    annotations=pickle.load(file)
    file.close()
    
    ipa.Gibbs_sampler_bio(df, annotations, Bio, noits=1000, burn=None, delta_bio=.1, all_out=False, zs=None)
    
    ex1=expected1[1]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[501]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[4]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[999]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
      

def test_Gibbs_sampler_bio_add():
    file = open('test_Gibbs_sampler_bio.pkl', 'rb')
    df=pickle.load(file)
    Bio=pickle.load(file)
    annotations=pickle.load(file)
    file.close()
    
    ipa.Gibbs_sampler_bio_add(df, annotations, Bio, noits=1000, burn=None, delta_bio=.1, delta_add=0.1,
                              all_out=False, zs=None)
    
    ex1=expected1[1]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[501]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[4]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=expected1[999]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
  

def test_simpleIPA():
    file = open('test_simpleIPA.pkl','rb')
    df=pickle.load(file)
    DB=pickle.load(file)
    adducts=pickle.load(file)
    dfMS2=pickle.load(file)
    DBMS2=pickle.load(file)
    Bio=pickle.load(file)
    expected1=pickle.load(file)
    expected2=pickle.load(file)
    expected3=pickle.load(file)
    expected4=pickle.load(file)
    expected5=pickle.load(file)
    file.close()
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=None, DBMS2=None, noits=1000,
                                burn=None, delta_add=None, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    assert(annotations[1].equals(expected1[1]))
    assert(annotations[501].equals(expected1[501]))
    assert(annotations[4].equals(expected1[4]))
    assert(annotations[999].equals(expected1[999]))
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=None, DBMS2=None, noits=1000,
                                burn=None, delta_add=None, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=2, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    assert(annotations[1].equals(expected1[1]))
    assert(annotations[501].equals(expected1[501]))
    assert(annotations[4].equals(expected1[4]))
    assert(annotations[999].equals(expected1[999]))

    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=None, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    assert(annotations[1].equals(expected2[1]))
    assert(annotations[501].equals(expected2[501]))
    assert(annotations[4].equals(expected2[4]))
    assert(annotations[999].equals(expected2[999]))
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=None, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=2, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    assert(annotations[1].equals(expected2[1]))
    assert(annotations[501].equals(expected2[501]))
    assert(annotations[4].equals(expected2[4]))
    assert(annotations[999].equals(expected2[999]))
    
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=expected3[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=3, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=expected3[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected3[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    

    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=None, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])
    
    ex3=expected4[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected4[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=None, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=2, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])
    
    ex3=expected4[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    

    ex3=expected4[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
 
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=expected5[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    
    ex3=expected5[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    
    annotations = ipa.simpleIPA(df, ionisation=1, DB=DB, adductsAll=adducts,
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=3, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=expected5[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=expected5[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))


def test_diff():
    assert(MS2compare.diff([20,10])==[-10])
