import pandas
import unicodedata

__author__ = "Francesco Del Carratore"
__maintainer__ = "Francesco Del Carratore"
__email__ = "francescodc87@gmail.com"

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




def get_annotations(id,IPA_data):
    if id not in IPA_data:   
        return None
    out = {}
    data = IPA_data[id]  
    out["IPA_id"] = ', '.join(str(i) for i in list(data["id"]))
    out["IPA_name"] = ', '.join(str(i) for i in list(data["name"]))
    out["IPA_formula"] = ', '.join(str(i) for i in list(data["formula"]))
    out["IPA_adduct"] = ', '.join(str(i) for i in list(data["adduct"]))
    out["IPA_mz"] = ', '.join(str(i) for i in list(data["m/z"]))
    out["IPA_charge"] = ', '.join(str(i) for i in list(data["charge"]))
    out["IPA_ppm"] = ', '.join(str(i) for i in list(data["ppm"]))
    out["IPA_isotope_pattern_score"] = ', '.join(str(i) for i in list(data["isotope pattern score"]))
    out["IPA_fragmentation_pattern_score"] = ', '.join(str(i) for i in list(data["fragmentation pattern score"]))
    out["IPA_prior"] = ', '.join(str(i) for i in list(data["prior"]))
    out["IPA_post"] = ', '.join(str(i) for i in list(data["post"]))
    out["IPA_post_Gibbs"] = ', '.join(str(i) for i in list(data["post Gibbs"]))
    out["IPA_post_chi_square_pval"] = ', '.join(str(i) for i in list(data["chi-square pval"]))
    
    check = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ,()1234567890-[]_.+-'"
    new_string = ""
    for i in out["IPA_name"]:
        if i in check:
            new_string += i
        else:
            try:
                sc = unicodedata.name(i)
                new_string += f"[{sc}]"
            except:
                print(f"{i} does not exist in unicode")
    out["IPA_name"] = new_string
    return out

def input_string(data,exept):
    if data == None:
        return None
    out = ""
    for i in data:
        if i not in exept:
            out += f'				<annotation unit="">\n					<label>{i}</label>\n					<value>{data[i]}</value>\n					<valuetype>STRING</valuetype>\n				</annotation>\n'
    return out