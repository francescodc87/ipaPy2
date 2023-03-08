from ipaPy2 import ipa
import pickle

def test_clusterFeatures():
    file = open('tests/test_clusterFeatures.pkl', 'rb')
    df= pickle.load(file)
    expected = pickle.load(file)
    file.close()
    df2 = ipa.clusterFeatures(df)
    assert(df2.equals(expected))
    