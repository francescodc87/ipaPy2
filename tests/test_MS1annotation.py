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