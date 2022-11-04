import pandas

__author__ = "Francesco Del Carratore"
__maintainer__ = "Francesco Del Carratore"
__email__ = "francesco.delcarratore@gmail.com"

def which(self):
    """equivalent of the which function in R 
    """
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)

def whichNone(self):
    """give the indices of which element is equal to None 
    """
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if x is None]
    return(indices)


def check_ded(f1,f2):
    """check if formula 2 is included in formula1 
    """
    A=pandas.DataFrame(f1.composition(), columns=['Atoms','num','a','b']).iloc[:,range(0,2)]
    B=pandas.DataFrame(f2.composition(), columns=['Atoms','num','a','b']).iloc[:,range(0,2)]
    AB=pandas.merge(A,B, on='Atoms', how='outer')
    AB=AB.fillna(0)
    check=AB.iloc[:,1]-AB.iloc[:,2]
    res=all(check>0)
    return(res)
