from ipaPy2 import MS2compare

def test_diff():
    assert(MS2compare.diff([20,10])==[-10])
    