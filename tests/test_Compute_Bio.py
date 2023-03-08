from ipaPy2 import ipa
import pickle

def test_Compute_Bio():
    file = open('tests/test_Compute_Bio.pkl', 'rb')
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
    
    
        