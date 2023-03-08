from ipaPy2 import ipa
import pickle

def test_simpleIPA():
    file = open('tests/test_simpleIPA.pkl','rb')
    df=pickle.load(file)
    DB=pickle.load(file)
    adducts=pickle.load(file)
    dfMS2=pickle.load(file)
    DBMS2=pickle.load(file)
    Bio=pickle.load(file)
    expected1=pickle.load(file)
    expected2=pickle.load(file)
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

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
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

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
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
    
    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
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
    
    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
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

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
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

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
    
    ex3=annotations[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown']))
