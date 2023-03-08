from ipaPy2 import ipa
import pickle


def test_Gibbs_sampler_add():
    file = open('tests/test_Gibbs_sampler_add.pkl', 'rb')
    df=pickle.load(file)
    annotations=pickle.load(file)
    file.close()
    
    ipa.Gibbs_sampler_add(df, annotations, noits=1000, burn=None, delta_add=.1, all_out=False, zs=None)
    
    ex1=annotations[1]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=annotations[501]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00148'])> \
    float(ex1['post Gibbs'][ex1['id']=='C00763'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=annotations[4]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    
    ex1=annotations[999]
    assert(float(ex1['post Gibbs'][ex1['id']=='C00079'])> \
    float(ex1['post Gibbs'][ex1['id']=='C02265'])> \
    float(ex1['post Gibbs'][ex1['id']=='Unknown']))
    