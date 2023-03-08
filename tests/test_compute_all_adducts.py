from ipaPy2 import ipa
import pickle

def test_compute_all_adducts():
    file = open('tests/test_compute_all_adducts.pkl', 'rb')
    adducts=pickle.load(file)
    DB=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    
    out1 = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=1)
    out2 = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=2)
    out3 = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=1)
    out4 = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=2)
     
    assert(out1.equals(expected1) and out2.equals(expected1) and out3.equals(expected2) and out4.equals(expected2))
