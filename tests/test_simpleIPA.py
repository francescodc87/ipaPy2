from ipaPy2 import ipa
import pickle

def test_simpleIPA1():
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

    assert(list(annotations[1]['post'])==list(expected1[1]['post']) and \
        list(annotations[1]['id'])==list(expected1[1]['id']))

    assert(list(annotations[501]['post'])==list(expected1[501]['post']) and \
        list(annotations[501]['id'])==list(expected1[501]['id']))

    assert(list(annotations[4]['post'])==list(expected1[4]['post']) and \
        list(annotations[4]['id'])==list(expected1[4]['id']))

    assert(list(annotations[999]['post'])==list(expected1[999]['post']) and \
        list(annotations[999]['id'])==list(expected1[999]['id']))


def test_simpleIPA2():
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
                            ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                            burn=None, delta_add=None, delta_bio=None,
                            Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                            ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                            ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                            pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                            connections=['C4H2'])

    assert(list(annotations[1]['post'])==list(expected2[1]['post']) and \
        list(annotations[1]['id'])==list(expected2[1]['id']))

    assert(list(annotations[501]['post'])==list(expected2[501]['post']) and \
        list(annotations[501]['id'])==list(expected2[501]['id']))

    assert(list(annotations[4]['post'])==list(expected2[4]['post']) and \
        list(annotations[4]['id'])==list(expected2[4]['id']))

    assert(list(annotations[999]['post'])==list(expected2[999]['post']) and \
        list(annotations[999]['id'])==list(expected2[999]['id']))

    
def test_simpleIPA3():
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
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=None,
                                Bio=None, mode='reactions', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    

def test_simpleIPA4():
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
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=None, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])
    
    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))


   
def test_simpleIPA5():
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
                                ppm=5, dfMS2=dfMS2, DBMS2=DBMS2, noits=1000,
                                burn=None, delta_add=0.1, delta_bio=0.1,
                                Bio=Bio, mode='connections', CSunk=0.5, isodiff=1,
                                ppmiso=100, ncores=1, me=0.000548579909065, ratiosd=0.9,
                                ppmunk=None, ratiounk=None, ppmthr=None, pRTNone=None,
                                pRTout=None, mzdCS=0, ppmCS=10, evfilt=False,
                                connections=['C4H2'])

    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
