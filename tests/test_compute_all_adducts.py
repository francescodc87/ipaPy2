from ipaPy2 import ipa
import pickle

def test_compute_all_adducts():
    file = open('tests/test_compute_all_adducts.pkl', 'rb')
    adducts=pickle.load(file)
    DB=pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    
    out = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=1)
    assert(out.equals(expected1))

    out = ipa.compute_all_adducts(adducts, DB, ionisation=1, ncores=2)
    assert(out.equals(expected1))
    
    out = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=1)
    assert(out.equals(expected2))

    out = ipa.compute_all_adducts(adducts, DB, ionisation=-1, ncores=2)
    assert(out.equals(expected2))