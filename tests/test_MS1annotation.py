from ipaPy2 import ipa
import pickle

def test_MS1annotation():
    file = open('tests/test_MS1annotation.pkl', 'rb')
    df=pickle.load(file)
    allAdds=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    annotations=ipa.MS1annotation(df,allAdds,ppm=5,
                                  me = 5.48579909065e-04,ratiosd=0.9,
                                  ppmunk=None,ratiounk=None,ppmthr=None,
                                  pRTNone=None,pRTout=None,ncores=1)

    annotations2=ipa.MS1annotation(df,allAdds,ppm=5,
                                  me = 5.48579909065e-04,ratiosd=0.9,
                                  ppmunk=None,ratiounk=None,ppmthr=None,
                                  pRTNone=None,pRTout=None,ncores=2)    
    
    assert(annotations[1].equals(expected1[1]) and annotations[501].equals(expected1[501]) and annotations[4].equals(expected1[4]) and annotations[999].equals(expected1[999]) and annotations2[1].equals(expected2[1]) and annotations2[501].equals(expected2[501]) and annotations2[4].equals(expected2[4]) and annotations2[999].equals(expected2[999]))