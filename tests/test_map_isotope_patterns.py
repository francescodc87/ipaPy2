from ipaPy2 import ipa
import pickle

def test_map_isotope_patterns():
    file = open('tests/test_map_isotope_patterns.pkl', 'rb')
    df= pickle.load(file)
    expected1 = pickle.load(file)
    expected2 = pickle.load(file)
    file.close()
    
    df_out =df.copy()
    ipa.map_isotope_patterns(df_out,ionisation=1)
    assert(df_out.equals(expected1))
    
    df_out =df.copy()
    ipa.map_isotope_patterns(df_out,ionisation=-1)
    assert(df_out.equals(expected2))