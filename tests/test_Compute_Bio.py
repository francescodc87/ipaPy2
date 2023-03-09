from ipaPy2 import ipa
import pickle


def test_Compute_Bio():
    file = open('tests/test_Compute_Bio.pkl', 'rb')
    DB=pickle.load(file)
    annotations=pickle.load(file)
    file.close()
    
    Bio1=ipa.Compute_Bio(DB,annotations,mode='reactions',ncores=1)
    assert(len(Bio1.index)==1 and sorted(list(Bio1.loc[0]==['C00079', 'C02265'])))