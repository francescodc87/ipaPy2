from ipaPy2 import ipa
import pickle

def test_MS1annotation():
    file = open('tests/test_MS1annotation.pkl', 'rb')
    df=pickle.load(file)
    allAdds=pickle.load(file)
    expected1 = pickle.load(file)
    file.close()
    
    annotations=ipa.MS1annotation(df,allAdds,ppm=5,
                                  me = 5.48579909065e-04,ratiosd=0.9,
                                  ppmunk=None,ratiounk=None,ppmthr=None,
                                  pRTNone=None,pRTout=None,ncores=1)

    assert(list(annotations[1]['post'])==list(expected1[1]['post']) and \
        list(annotations[1]['id'])==list(expected1[1]['id']))

    assert(list(annotations[501]['post'])==list(expected1[501]['post']) and \
        list(annotations[501]['id'])==list(expected1[501]['id']))

    assert(list(annotations[4]['post'])==list(expected1[4]['post']) and \
        list(annotations[4]['id'])==list(expected1[4]['id']))

    assert(list(annotations[999]['post'])==list(expected1[999]['post']) and \
        list(annotations[999]['id'])==list(expected1[999]['id']))