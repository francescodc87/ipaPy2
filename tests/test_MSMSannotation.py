from ipaPy2 import ipa
import pickle

def test_MSMSannotation():
    file = open('tests/test_MS2annotation.pkl', 'rb')
    df=pickle.load(file)
    allAdds=pickle.load(file)
    dfMS2=pickle.load(file)
    DBMS2=pickle.load(file)
    expected1 = pickle.load(file)
    file.close()
    
    annotations=ipa.MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm=5,
                               me=0.000548579909065,ratiosd=0.9,
                               ppmunk=None,ratiounk=None,ppmthr=None,
                               pRTNone=None,pRTout=None,mzdCS=0,ppmCS=10,
                               CSunk=0.7,evfilt=False,ncores=1)

    assert(list(annotations[1]['post'])==list(expected1[1]['post']) and \
        list(annotations[1]['id'])==list(expected1[1]['id']))

    assert(list(annotations[501]['post'])==list(expected1[501]['post']) and \
        list(annotations[501]['id'])==list(expected1[501]['id']))

    assert(list(annotations[4]['post'])==list(expected1[4]['post']) and \
        list(annotations[4]['id'])==list(expected1[4]['id']))

    assert(list(annotations[999]['post'])==list(expected1[999]['post']) and \
        list(annotations[999]['id'])==list(expected1[999]['id']))