from ipaPy2 import ipa
import pickle


def test_Gibbs_sampler_add():
    file = open('tests/test_Gibbs_sampler_add.pkl', 'rb')
    df=pickle.load(file)
    annotations=pickle.load(file)
    file.close()
    
    ipa.Gibbs_sampler_add(df, annotations, noits=1000, burn=None, delta_add=.1, all_out=False, zs=None)
    
    ex3=annotations[1]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[501]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00148'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C00763'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[4]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))
    
    ex3=annotations[999]
    assert(float(ex3['post Gibbs'][ex3['id']=='C00079'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='C02265'].iloc[0])> \
    float(ex3['post Gibbs'][ex3['id']=='Unknown'].iloc[0]))