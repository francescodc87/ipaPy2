from ipaPy2 import ipa
import pickle


def test_Compute_Bio():
    file = open('tests/test_Compute_Bio.pkl', 'rb')
    DB=pickle.load(file)
    annotations=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    Bio1=ipa.Compute_Bio(DB,annotations,mode='reactions',ncores=1)
    Bio2=ipa.Compute_Bio(DB,annotations,mode='reactions',ncores=2)
    Bio3=ipa.Compute_Bio(DB,annotations,mode='connections',connections=['C4H2'],ncores=1)
    Bio4=ipa.Compute_Bio(DB,annotations,mode='connections',connections=['C4H2'],ncores=2)
    assert(Bio1.equals(expected1) and Bio2.equals(expected1) and Bio3.equals(expected2) and Bio4.equals(expected2))