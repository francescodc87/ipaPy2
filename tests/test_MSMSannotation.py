from ipaPy2 import ipa
import pickle

def test_MSMSannotation():
    file = open('tests/test_MS2annotation.pkl', 'rb')
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
    annotations2=ipa.MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm=5,
                               me=0.000548579909065,ratiosd=0.9,
                               ppmunk=None,ratiounk=None,ppmthr=None,
                               pRTNone=None,pRTout=None,mzdCS=0,ppmCS=10,
                               CSunk=0.7,evfilt=False,ncores=2)

    assert(annotations[1].equals(expected1[1]) and annotations[501].equals(expected1[501]) and annotations[4].equals(expected1[4]) and annotations[999].equals(expected1[999]) and annotations[1].equals(expected1[1]) and annotations[501].equals(expected1[501]) and annotations[4].equals(expected1[4]) and annotations[999].equals(expected1[999]))