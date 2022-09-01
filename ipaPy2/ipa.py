from __future__ import annotations
import pandas
import numpy
import time
import molmass
from scipy import stats
import math
import random
import collections
import itertools
import multiprocessing
from functools import partial
import statistics
from tqdm import tqdm
from ipaPy2 import MS2compare
from ipaPy2 import util






def map_isotope_patterns(df,isoDiff=1, ppm=100, ionisation=1):
    """mapping isotope patterns in MS1 data.
    Inputs:
        df: pandas dataframe with the following columns:
            1) ids: an unique id for each feature
            2) rel.ids: relation ids. In a previous step of the data processing pipeline,features are
                        clustered based on peak shape similarity/retention time. Features in the same
                        cluster are likely to come from the same metabolite. All isotope patterns must
                        be in the same rel.id cluster.
            3) mzs: mass-to-charge ratios, usually the average across different samples.
            4) RTs: retention times in seconds, usually the average across different samples.
            5) maxInts: maximum intensity detected for each feature across samples (either peak area
                        or peak intensity)
        isoDiff: Default value 1. Difference between isotopes of charge 1, does not need to be exact
        ppm: Default value 100. Maximum ppm value allowed between 2 isotopes. It is very high on
             purpose
        ionisation: Default value 1. positive = 1, negative = -1
    Output:
         df: the main input is modified by adding and populating the following colums
             6) relationship: the possible values are:
                              - bp: basepeak, most intense peak within each realtion id
                              - bp|isotope: isotope of the basepeak
                              - potential bp: most intense peak within each isotope pattern
                                              (excluding the basepeak)
                              - potential bp|isotope: isotope of one potential bp
             7) isotope pattern: feature used to cluster the different isotope patterns within
                                 the same relation id
             8) charge: predicted charge based on the isotope pattern (1,2,3 or -1,-2,-3 are
                        the only values considered)
    """
    print("mapping isotope patterns ....")
    start = time.time()
    if isinstance(df, pandas.DataFrame):
        relIds = df.iloc[:,1]
        relIds = list(set(relIds))
        df['relationship'] = [None] * len(df.index)
        df['isotope pattern'] = [None] * len(df.index)
        df['charge'] = [None] * len(df.index)
        df['ind'] = list(range(0,len(df.index)))
        f1 = False
        f2 = False
        for g in relIds:
            ind = util.which(df.iloc[:,1] == g)
            dfg = df.iloc[ind,:].copy()
            dfg = dfg.sort_values(by=['mzs'])
            c = 0
            f1=False
            f2=False
            for k in range(0,len(dfg.index)-1):
                ppm1 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+isoDiff))/(dfg.iloc[k,2]+isoDiff))*(10**6))
                ppm2 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/2)))/(dfg.iloc[k,2]+(isoDiff/2)))*(10**6))
                ppm3 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/3)))/(dfg.iloc[k,2]+(isoDiff/3)))*(10**6))
                ppm4 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/4)))/(dfg.iloc[k,2]+(isoDiff/3)))*(10**6))
                ppm5 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/5)))/(dfg.iloc[k,2]+(isoDiff/3)))*(10**6))
                indiso1 = util.which(ppm1 <= ppm)
                indiso2 = util.which(ppm2 <= ppm)
                indiso3 = util.which(ppm3 <= ppm)
                indiso4 = util.which(ppm4 <= ppm)
                indiso5 = util.which(ppm5 <= ppm)
                if len(indiso5) >0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso5,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso5,6] = c
                    dfg.iloc[k,7] = 3*ionisation
                    dfg.iloc[indiso5,7] = 3*ionisation
                    f2=True
                elif len(indiso4) > 0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso4,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso4,6] = c
                    dfg.iloc[k,7] = 3*ionisation
                    dfg.iloc[indiso4,7] = 3*ionisation
                    f2=True
                elif len(indiso3) > 0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso3,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso3,6] = c
                    dfg.iloc[k,7] = 3*ionisation
                    dfg.iloc[indiso3,7] = 3*ionisation
                    f2=True
                elif len(indiso2)>0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso2,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso2,6] = c
                    dfg.iloc[k,7] = 2*ionisation
                    dfg.iloc[indiso2,7] = 2*ionisation
                    f2=True
                elif len(indiso1)>0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso1,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso1,6] = c
                    dfg.iloc[k,7] = 1*ionisation
                    dfg.iloc[indiso1,7] = 1*ionisation
                    f2=True
                else:
                    f1 = dfg.iloc[k,5]!=None
                if f1 and f2:
                    c=c+1
                    f1=False
                    f2=False

            df.iloc[dfg.iloc[:,8],6] = dfg.iloc[:,6]
            df.iloc[dfg.iloc[:,8],5] = dfg.iloc[:,5]
            df.iloc[dfg.iloc[:,8],7] = dfg.iloc[:,7]

            ##fix isotope patterns
            iso_patterns = list(set(dfg.iloc[:,6]))
            for k in iso_patterns:
                if k!=None:
                    indpbpiso =dfg.iloc[util.which(dfg.iloc[:,6]==k),8]
                    tmp = df.iloc[indpbpiso,:].copy()
                    indpbp= tmp.iloc[util.which(tmp.iloc[:,4]==max(tmp.iloc[:,4]))[0],8]

                    df.iloc[indpbpiso,5] = "potential bp|isotope"
                    df.iloc[indpbp,5] = "potential bp"



            #find basepeaks
            indbp=dfg.iloc[util.which(dfg.iloc[:,4]==max(dfg.iloc[:,4]))[0],8]
            indbpiso =dfg.iloc[util.which(dfg.iloc[:,6]==df.iloc[indbp,6]),8]
            df.iloc[indbpiso,5] = "bp|isotope"
            df.iloc[indbp,5] = "bp"
        df.drop(columns=['ind'],inplace=True)
    else:
        raise Exception("""'map_isotope_patterns' method can only be applied to pandas dataframe.""")
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')



def compute_all_adducts(adducts_file, DB):
    """compute all adducts table based on the information present in the database
    Inputs:
        adducts_file: path to a csv file containing inforamtion on all possible adducts. The file
                      must be in the following format, and column names must be the same:
        Name    calc	    Charge	Mult	Mass	    Ion_mode	Formula_add	Formula_ded	Multi
        M+H     M+1.007276	1	    1	    1.007276	positive	H1	        FALSE	    1
        M-H	    M-1.007276	-1	    1	    -1.007276	negative	FALSE	    H1	        1
        DB: pandas dataframe containing the database against which the annotation is performed. The DB must
            contain the following columns in this exact order (optional fields can contain None):
                - id: unique id of the database entry (e.g., 'C00031') - necessary
                - name: compund name (e.g., 'D-Glucose') - necessary
                - formula: chemical formula (e.g., 'C6H12O6') - necessary
                - inchi: inchi string - optional
                - smiles: smiles string - optional
                - RT: if known, retention time range (in seconds) where this compound is expected
                      to elute (e.g., '30;60') - optional
                - adducts: list of adducts that should be considered for this entry (e.g.,'M+Na;M+H;M+')
                - description: comments on the entry - optional
                - pk: previous knowledge on the likelihood of this compoud to be present in the sample
                      analyse. The value has to be between 1 (compound likely to be present in the sample)
                      and 0 (compound cannot be present in the sample).
                - MS2: id for the MS2 database entries related to this compound (optional)
                - reactions: list of reactions ids involving this compound (e.g., 'R00010 R00015 R00028')-optional 
    Output:
         allAdds: pandas dataframe containing the information on all the possible adducts given the database.
                  Example:
        id	    name	    adduct	formula	    charge	m/z	        RT	    pk	    MS2
        C00031	D-Glucose	M+H	    C6H13O6	    1	    181.070664	30;60	1.0	    EMBL-MCF_spec393569_1
    """
    print("computing all adducts ....")
    start = time.time()
    adductsAll = pandas.read_csv(adducts_file)
    data=[]
    for db in range(0,len(DB.index)):
        data.append(all_adducts_iter(DB,adductsAll,db))
    allAdds =pandas.concat(data,ignore_index=True)
    allAdds.columns=['id','name','adduct','formula','charge','m/z','RT','pk','MS2']
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(allAdds)


def all_adducts_iter(DB,adductsAll,db):
    f = DB['formula'][db]
    f = molmass.Formula(f)
    M=f.isotope.mass
    recs = []
    adds = DB['adducts'][db]
    adds = adds.split(';')
    adducts = adductsAll[adductsAll['Name'].isin(adds)]
    for add in range(0,len(adducts.index)):
        f1=f
        id=DB['id'][db]
        name=DB['name'][db]
        addtype=adducts.iloc[add,0]
        charge = adducts.iloc[add,2]
        rt = DB['RT'][db]
        pk = DB['pk'][db]
        ms2ind = DB['MS2'][db]
        mz=M/abs(charge)
        mz=mz*adducts.iloc[add,3]
        mz=mz+adducts.iloc[add,4]
        if adducts.iloc[add,6] != 'FALSE':
            f1=f1.__mul__(int(adducts.iloc[add,8]))
            f1=f1.__add__(molmass.Formula(adducts.iloc[add,6]))
        if adducts.iloc[add,7] != 'FALSE':
            ###I need to add a check here to see if the sub is possible
            f1=f1.__mul__(int(adducts.iloc[add,8]))
            sub=molmass.Formula(adducts.iloc[add,7])
            if util.check_ded(f1,sub): 
                f1=f1.__sub__(sub)
        form=f1.formula    
        recs.append([id,name,addtype,form,charge,mz,rt,pk,ms2ind])
    return(pandas.DataFrame(recs)) 
    
    

def compute_all_adducts_Parallel(adducts_file, DB, ncores=1):
    """compute all adducts table based on the information present in the database - parallelized version
    Inputs:
        adducts_file: path to a csv file containing inforamtion on all possible adducts. The file
                      must be in the following format, and column names must be the same:
        Name    calc	    Charge	Mult	Mass	    Ion_mode	Formula_add	Formula_ded	Multi
        M+H     M+1.007276	1	    1	    1.007276	positive	H1	        FALSE	    1
        M-H	    M-1.007276	-1	    1	    -1.007276	negative	FALSE	    H1	        1
        DB: pandas dataframe containing the database against which the annotation is performed. The DB must
            contain the following columns in this exact order (optional fields can contain None):
                - id: unique id of the database entry (e.g., 'C00031') - necessary
                - name: compund name (e.g., 'D-Glucose') - necessary
                - formula: chemical formula (e.g., 'C6H12O6') - necessary
                - inchi: inchi string - optional
                - smiles: smiles string - optional
                - RT: if known, retention time range (in seconds) where this compound is expected
                      to elute (e.g., '30;60') - optional
                - adducts: list of adducts that should be considered for this entry (e.g.,'M+Na;M+H;M+')
                - description: comments on the entry - optional
                - pk: previous knowledge on the likelihood of this compoud to be present in the sample
                      analyse. The value has to be between 1 (compound likely to be present in the sample)
                      and 0 (compound cannot be present in the sample).
                - MS2: id for the MS2 database entries related to this compound (optional)
                - reactions: list of reactions ids involving this compound (e.g., 'R00010 R00015 R00028')-optional
        ncores: default value 1. Number of cores used
    Output:
         allAdds: pandas dataframe containing the information on all the possible adducts given the database.
                  Example:
        id	    name	    adduct	formula	    charge	m/z	        RT	    pk	    MS2
        C00031	D-Glucose	M+H	    C6H13O6	    1	    181.070664	30;60	1.0	    EMBL-MCF_spec393569_1
    """
    print("computing all adducts ....")
    start = time.time()
    adductsAll = pandas.read_csv(adducts_file)
    pool_obj = multiprocessing.Pool(ncores)
    data = pool_obj.map(partial(all_adducts_iter,DB,adductsAll),range(0,len(DB.index)))
    pool_obj.terminate()
    allAdds =pandas.concat(data,ignore_index=True)
    allAdds.columns=['id','name','adduct','formula','charge','m/z','RT','pk','MS2']
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(allAdds)


def MS1annotation(df,allAdds,ppm,me = 5.48579909065e-04,ratiosd=0.9,ppmunk=None,ratiounk=None,ppmthr=None, pRTNone=None, pRTout=None):
    """Annotation of the dataset base on the MS1 information. Prior probabilities are based on mass only, while post probabilities
       are based on mass, RT, previous knowledge and isotope patterns.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        allAdds: pandas dataframe containing the information on all the possible adducts given the database. It should
                 be the output of either ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
        ppm: accuracy of the MS instrument used
        me: accurate mass of the electron. Default 5.48579909065e-04
        ratiosd: default 0.9. It represents the acceptable ratio between predicted intensity and observed intesity of isotopes.
                 it is used to compute the shape parameters of the lognormal distribution used to calculate the isotope pattern
                 scores as sqrt(1/ratiosd)
        ppmunk: ppm associated to the 'unknown' annotation. If not provided equal to ppm.
        ratiounk: isotope ratio associated to the 'unknown' annotation. If not provided equal to 0.5
        ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to 2*ppm
        pRTNone: Multiplicative factor for the RT if no RTrange present in the database. If not provided equal to 0.8
        pRTout: Multiplicative factor for the RT if measured RT is outside the RTrange present in the database.
                If not provided equal to 0.4
    Output:
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
    """
    print("annotating based on MS1 information....")
    start = time.time()
    if ppmunk is None:
        ppmunk = ppm
    if ppmthr is None:
        ppmthr = 2*ppm
    if ratiounk is None:
        ratiounk = 0.5
    if pRTNone is None:
        pRTNone = 0.8
    if pRTout is None:
        pRTout = 0.4
    annotations={}
    ind = util.which(df.iloc[:,5]=='bp')+ util.which(df.iloc[:,5]=='potential bp') + util.whichNone(df.iloc[:,5])
    ind.sort()
    sigmaln = math.sqrt(1/ratiosd)
    data=[]
    for k in ind:
        data.append(MS1_ann_iter(df,allAdds,ppm,me,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,sigmaln,k))
    keys = list(df.iloc[ind,0])
    annotations = dict(zip(keys, data))
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(annotations)        


def MS1annotation_Parallel(df,allAdds,ppm,me = 5.48579909065e-04,ratiosd=0.9,ppmunk=None, ratiounk=None,ppmthr=None, pRTNone=None, pRTout=None, ncores=1):
    """Annotation of the dataset base on the MS1 information. Prior probabilities are based on mass only, while post probabilities
       are based on mass, RT, previous knowledge and isotope patterns - parallelized version.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        allAdds: pandas dataframe containing the information on all the possible adducts given the database. It should
                 be the output of either ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
        ppm: accuracy of the MS instrument used
        me: accurate mass of the electron. Default 5.48579909065e-04
        ratiosd: default 0.9. It represents the acceptable ratio between predicted intensity and observed intesity of isotopes.
                 it is used to compute the shape parameters of the lognormal distribution used to calculate the isotope pattern
                 scores as sqrt(1/ratiosd)
        ppmunk: ppm associated to the 'unknown' annotation. If not provided equal to ppm.
        ratiounk: isotope ratio associated to the 'unknown' annotation. If not provided equal to 0.5
        ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to 2*ppm
        pRTNone: Multiplicative factor for the RT if no RTrange present in the database. If not provided equal to 0.8
        pRTout: Multiplicative factor for the RT if measured RT is outside the RTrange present in the database.
                If not provided equal to 0.4
        ncores: default value 1. Number of cores used
    Output:
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
    """
    print("annotating based on MS1 information....")
    print("annotating and computing isotope pattern scores ....")
    start = time.time()
    if ppmunk is None:
        ppmunk = ppm
    if ppmthr is None:
        ppmthr = 2*ppm
    if ratiounk is None:
        ratiounk = 0.5
    if pRTNone is None:
        pRTNone = 0.8
    if pRTout is None:
        pRTout = 0.4
    
    annotations={}
    ind = util.which(df.iloc[:,5]=='bp')+ util.which(df.iloc[:,5]=='potential bp') + util.whichNone(df.iloc[:,5])
    ind.sort()
    sigmaln = math.sqrt(1/ratiosd)
    pool_obj = multiprocessing.Pool(ncores)
    data = pool_obj.map(partial(MS1_ann_iter,df,allAdds,ppm,me,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,sigmaln),ind)
    pool_obj.terminate()
    keys = list(df.iloc[ind,0])
    annotations = dict(zip(keys, data))
    ##convert data into dictionary! annotations[df.iloc[k,0]]=tmp
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(annotations)        


def MS1_ann_iter(df,allAdds,ppm,me,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,sigmaln,k):
    mm= df.iloc[k,2]
    rtm=df.iloc[k,3]
    ppms = ((mm-allAdds.iloc[:,5])/allAdds.iloc[:,5])*(10**6) 
    hits = util.which(abs(ppms)<=ppmthr)
    priors = list(ppms[hits].copy())
    priors.append(ppmunk)
    priors = stats.norm(0, ppm/2).pdf(priors)
    priors=priors/sum(priors)
    pks = list(allAdds.pk[hits].copy())
    pks.append(1)
    if len(hits)>0: #have to add the unknowm
        tmp = allAdds.iloc[hits,[0,1,3,2,5,4,6]].copy()
        tmp = tmp.rename(columns={'RT':'RT range'})
        tmp['ppm'] = ppms[hits]
        tmp['isotope pattern score'] = [None]*len(hits)
        tmp['fragmentation pattern score'] = [None]*len(hits)
        tmp['prior'] = [None]*len(hits)
        tmp['post'] = [None]*len(hits)
        uk = pandas.DataFrame({'id':["Unknown"],'name':["Unknown"],'formula':[None],'adduct':[None],'m/z':[None],
                       'charge':[None],'RT range':[None],'ppm':[ppm],
                       'isotope pattern score':[None],
                       'fragmentation pattern score':[None],'prior':[None],'post':[None]})
        tmp=pandas.concat([tmp,uk])
        ### compute isotope pattern scores
        relid = list(df.iloc[util.which(df.iloc[:,0]==df.iloc[k,0]),1])[0]
        isoid = list(df.iloc[util.which(df.iloc[:,0]==df.iloc[k,0]),6])[0]
        tmp.iloc[:,10] = priors
        if isoid is not None:
            indiso = list(set(util.which(df.iloc[:,1]==relid)) & set(util.which(df.iloc[:,6]==isoid)))
            ISm = df.iloc[indiso,[2,4]]
            ### computing isotope patterns scores
            pisoM=[]
            pisoI=[]
            for isp in range(0,len(hits)):
                mtmp = molmass.Formula(tmp.iloc[isp,2])
                ch = tmp.iloc[isp,5]
                ISt = mtmp.spectrum()
                ISt = pandas.DataFrame(ISt).transpose()
                if(len(ISt.index)<len(ISm.index)):
                    pisoM.append(0)
                    pisoI.append(0)
                else:
                    if ch>0:
                        ISt.iloc[:,0] = (ISt.iloc[:,0]/abs(ch)) - me
                    else:
                        ISt.iloc[:,0] = (ISt.iloc[:,0]/abs(ch)) + me

                    ISt=ISt.sort_values(by=[0])
                    ISm=ISm.sort_values(by=['mzs'])
                    ISm.iloc[:,1] = ISm.iloc[:,1]/max(ISm.iloc[:,1])
                    ISt = ISt.iloc[range(0,len(ISm.index)),:]
                    ISt.iloc[:,1] = ISt.iloc[:,1]/max(ISt.iloc[:,1])
                    pMs = []
                    pIs = []
                    for ps in range(1,len(ISm.index)):
                        ppmk =((ISt.iloc[ps,0]-ISm.iloc[ps,0])/ISt.iloc[ps,0])*(10**6)
                        pMs.append(stats.norm(0, ppm/2).pdf(ppmk)/stats.norm(0, ppm/2).pdf(0))
                        pIs.append(stats.lognorm(scale=ISt.iloc[ps,1],s=sigmaln).pdf(ISm.iloc[ps,1])/stats.lognorm(scale=ISt.iloc[ps,1],s=sigmaln).pdf(math.exp(ISt.iloc[ps,1]-(sigmaln)**2) ))

                    pisoM.append(sum(pMs))
                    pisoI.append(sum(pIs))
            ###adding scores for unknown
            pisoM.append(sum([stats.norm(0, ppm).pdf(2*ppm)/stats.norm(0, ppm).pdf(0)]*len(ISm.index)))
            pisoI.append(sum([stats.lognorm(scale=ISt.iloc[0,1],s=sigmaln).pdf(ISt.iloc[0,1]+ratiounk*ISt.iloc[0,1])/stats.lognorm(scale=ISt.iloc[0,1],s=sigmaln).pdf(math.exp(ISt.iloc[0,1]-(sigmaln)**2))]*len(ISm.index)))    
            pisoM = pisoM/sum(pisoM)
            pisoI = pisoI/sum(pisoI)
            piso = pisoM*pisoI
            tmp.iloc[:,8] = piso/sum(piso)
            
        ### computing RT score
        pRT = [pRTNone]*len(tmp.index)
        RTranges = list(tmp['RT range'])

        for r in range(0,len(RTranges)):
            if RTranges[r]!=None:
                rtrange = RTranges[r].split(';')
                rtrange = [float(i) for i in rtrange]
                if rtm >= rtrange[0] and rtm <= rtrange[1]:
                    pRT[r]=1
                else:
                    pRT[r]=pRTout


        ### computing posteriors integrating iso score and pk and RTscore
        p1 = tmp.iloc[:,8] 
        p1 = [1 if v is None else v for v in p1]# isopattern score
        pks = [1 if v is None else v for v in pks]# previous knowledge score
        post = [a*b*c*d for a,b,c,d in zip(priors,p1,pks,pRT)]
        post = [x / sum(post) for x in post]
        tmp.iloc[:,11] = post
        tmp = tmp.sort_values(by=['post'], ascending=False)
    else:
        tmp=pandas.DataFrame({'id':["Unknown"],'name':["Unknown"],'formula':[None],'adduct':[None],'m/z':[None],
                       'charge':[None],'RT range':[None],'ppm':[ppm],'isotope pattern score':[None],
                       'fragmentation pattern score':[None],'prior':[1],'post':[1]})
    tmp.index=range(0,len(tmp.index))    
    return(tmp)


def MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm,me = 5.48579909065e-04,ratiosd=0.9,ppmunk=None, ratiounk=None,
ppmthr=None,pRTNone=None, pRTout=None,mzdCS=0, ppmCS=10, CSunk=0.5):
    """Annotation of the dataset base on the MS1 and MS2 information. Prior probabilities are based on mass only, while post probabilities
       are based on mass, RT, previous knowledge and isotope patterns.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        dfMS2: pandas dataframe containing the MS2 data. It must contain 3 columns
                -id: an unique id for each feature for which the MS2 spectrum was acquired (same as in df)
                -spectrum: string containing the spectrum inforamtion in the following format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'
                -ev: collision energy used to aquire the fragmentation spectrum
        allAdds: pandas dataframe containing the information on all the possible adducts given the database. It should
                 be the output of either ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
        DBMS2: pandas dataframe containing the database containing the MS2 information
        ppm: accuracy of the MS instrument used
        me: accurate mass of the electron. Default 5.48579909065e-04
        ratiosd: default 0.9. It represents the acceptable ratio between predicted intensity and observed intesity of isotopes.
                 it is used to compute the shape parameters of the lognormal distribution used to calculate the isotope pattern
                 scores as sqrt(1/ratiosd)
        ppmunk: ppm associated to the 'unknown' annotation. If not provided equal to ppm.
        ratiounk: isotope ratio associated to the 'unknown' annotation. If not provided equal to 0.5
        ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to 2*ppm
        pRTNone: Multiplicative factor for the RT if no RTrange present in the database. If not provided equal to 0.8
        pRTout: Multiplicative factor for the RT if measured RT is outside the RTrange present in the database.
                If not provided equal to 0.4
        mzdCS: maximum mz difference allowed when computing cosine similarity scores. If one wants to use this parameter
               instead of ppmCS, this must be set to 0. Default 0.
        ppmCS: maximum ppm allowed when computing cosine similarity scores. If one wants to use this parameter
               instead of mzdCS, this must be set to 0. Default 10.
        CSunk: cosine similarity score associated with the 'unknown' annotation. Default 0.3
    Output:
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
    """
    print("annotating based on MS1 and MS2 information....")
    start = time.time()
    if ppmunk is None:
        ppmunk = ppm
    if ppmthr is None:
        ppmthr = 2*ppm
    if ratiounk is None:
        ratiounk = 0.5
    if pRTNone is None:
        pRTNone = 0.8
    if pRTout is None:
        pRTout = 0.4
    annotations={}
    ind = util.which(df.iloc[:,5]=='bp')+ util.which(df.iloc[:,5]=='potential bp') + util.whichNone(df.iloc[:,5])
    ind.sort()
    sigmaln = math.sqrt(1/ratiosd)
    data=[]
    for k in ind:
        data.append(MSMS_ann_iter(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k))
    keys = list(df.iloc[ind,0])
    annotations = dict(zip(keys, data))
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(annotations)       


def MSMSannotation_Parallel(df,dfMS2,allAdds,DBMS2,ppm,me = 5.48579909065e-04,ratiosd=0.9,ppmunk=None, ratiounk=None,
ppmthr=None,pRTNone=None, pRTout=None,mzdCS=0, ppmCS=10, CSunk=0.5,ncores=1):
    """Annotation of the dataset base on the MS1 and MS2 information. Prior probabilities are based on mass only, while post probabilities
       are based on mass, RT, previous knowledge and isotope patterns - parallelized version.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        dfMS2: pandas dataframe containing the MS2 data. It must contain 3 columns
                -id: an unique id for each feature for which the MS2 spectrum was acquired (same as in df)
                -spectrum: string containing the spectrum inforamtion in the following format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'
                -ev: collision energy used to aquire the fragmentation spectrum
        allAdds: pandas dataframe containing the information on all the possible adducts given the database. It should
                 be the output of either ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
        DBMS2: pandas dataframe containing the database containing the MS2 information
        ppm: accuracy of the MS instrument used
        me: accurate mass of the electron. Default 5.48579909065e-04
        ratiosd: default 0.9. It represents the acceptable ratio between predicted intensity and observed intesity of isotopes.
                 it is used to compute the shape parameters of the lognormal distribution used to calculate the isotope pattern
                 scores as sqrt(1/ratiosd)
        ppmunk: ppm associated to the 'unknown' annotation. If not provided equal to ppm.
        ratiounk: isotope ratio associated to the 'unknown' annotation. If not provided equal to 0.5
        ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to 2*ppm
        pRTNone: Multiplicative factor for the RT if no RTrange present in the database. If not provided equal to 0.8
        pRTout: Multiplicative factor for the RT if measured RT is outside the RTrange present in the database.
                If not provided equal to 0.4
        mzdCS: maximum mz difference allowed when computing cosine similarity scores. If one wants to use this parameter
               instead of ppmCS, this must be set to 0. Default 0.
        ppmCS: maximum ppm allowed when computing cosine similarity scores. If one wants to use this parameter
               instead of mzdCS, this must be set to 0. Default 10.
        CSunk: cosine similarity score associated with the 'unknown' annotation. Default 0.3
        ncores: default value 1. Number of cores used
    Output:
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
    """
    print("annotating based on MS1 and MS2 information....")
    start = time.time()
    if ppmunk is None:
        ppmunk = ppm
    if ppmthr is None:
        ppmthr = 2*ppm
    if ratiounk is None:
        ratiounk = 0.5
    if pRTNone is None:
        pRTNone = 0.8
    if pRTout is None:
        pRTout = 0.4
    annotations={}
    ind = util.which(df.iloc[:,5]=='bp')+ util.which(df.iloc[:,5]=='potential bp') + util.whichNone(df.iloc[:,5])
    ind.sort()
    sigmaln = math.sqrt(1/ratiosd)
    pool_obj = multiprocessing.Pool(ncores)
    data = pool_obj.map(partial(MSMS_ann_iter,df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln),ind)
    pool_obj.terminate()
    keys = list(df.iloc[ind,0])
    annotations = dict(zip(keys, data))
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(annotations)        

def MSMS_ann_iter(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k):
    mm= df.iloc[k,2]
    rtm=df.iloc[k,3]
    ppms = ((mm-allAdds.iloc[:,5])/allAdds.iloc[:,5])*(10**6) 
    hits = util.which(abs(ppms)<=ppmthr)
    priors = list(ppms[hits].copy())
    priors.append(ppmunk)
    priors = stats.norm(0, ppm/2).pdf(priors)
    pks = list(allAdds.pk[hits].copy())
    pks.append(1)
    priors=priors/sum(priors)
    if len(hits)>0: #have to add the unknowm
        tmp = allAdds.iloc[hits,[0,1,3,2,5,4,6]].copy()
        tmp = tmp.rename(columns={'RT':'RT range'})
        tmp['ppm'] = ppms[hits]
        tmp['isotope pattern score'] = [None]*len(hits)
        tmp['fragmentation pattern score'] = [None]*len(hits)
        tmp['prior'] = [None]*len(hits)
        tmp['post'] = [None]*len(hits)
        uk = pandas.DataFrame({'id':["Unknown"],'name':["Unknown"],'formula':[None],'adduct':[None],'m/z':[None],
                       'charge':[None],'RT range':[None],'ppm':[ppm],
                       'isotope pattern score':[None],
                       'fragmentation pattern score':[None],'prior':[None],'post':[None]})
        tmp=pandas.concat([tmp,uk])
        ### compute isotope pattern scores
        relid = list(df.iloc[util.which(df.iloc[:,0]==df.iloc[k,0]),1])[0]
        isoid = list(df.iloc[util.which(df.iloc[:,0]==df.iloc[k,0]),6])[0]
        tmp.iloc[:,10] = priors
        if isoid is not None:
            indiso = list(set(util.which(df.iloc[:,1]==relid)) & set(util.which(df.iloc[:,6]==isoid)))
            ISm = df.iloc[indiso,[2,4]]
            ### computing isotope patterns scores
            pisoM=[]
            pisoI=[]
            for isp in range(0,len(hits)):
                mtmp = molmass.Formula(tmp.iloc[isp,2])
                ch = tmp.iloc[isp,5]
                ISt = mtmp.spectrum()
                ISt = pandas.DataFrame(ISt).transpose()
                if(len(ISt.index)<len(ISm.index)):
                    pisoM.append(0)
                    pisoI.append(0)
                else:
                    if ch>0:
                        ISt.iloc[:,0] = (ISt.iloc[:,0]/abs(ch)) - me
                    else:
                        ISt.iloc[:,0] = (ISt.iloc[:,0]/abs(ch)) + me

                    ISt=ISt.sort_values(by=[0])
                    ISm=ISm.sort_values(by=['mzs'])
                    ISm.iloc[:,1] = ISm.iloc[:,1]/max(ISm.iloc[:,1])
                    ISt = ISt.iloc[range(0,len(ISm.index)),:]
                    ISt.iloc[:,1] = ISt.iloc[:,1]/max(ISt.iloc[:,1])
                    pMs = []
                    pIs = []
                    for ps in range(1,len(ISm.index)):
                        ppmk =((ISt.iloc[ps,0]-ISm.iloc[ps,0])/ISt.iloc[ps,0])*(10**6)
                        pMs.append(stats.norm(0, ppm/2).pdf(ppmk)/stats.norm(0, ppm/2).pdf(0))
                        pIs.append(stats.lognorm(scale=ISt.iloc[ps,1],s=sigmaln).pdf(ISm.iloc[ps,1])/stats.lognorm(scale=ISt.iloc[ps,1],s=sigmaln).pdf(math.exp(ISt.iloc[ps,1]-(sigmaln)**2) ))

                    pisoM.append(sum(pMs))
                    pisoI.append(sum(pIs))
            ###adding scores for unknown
            pisoM.append(sum([stats.norm(0, ppm).pdf(2*ppm)/stats.norm(0, ppm).pdf(0)]*len(ISm.index)))
            pisoI.append(sum([stats.lognorm(scale=ISt.iloc[0,1],s=sigmaln).pdf(ISt.iloc[0,1]+ratiounk*ISt.iloc[0,1])/stats.lognorm(scale=ISt.iloc[0,1],s=sigmaln).pdf(math.exp(ISt.iloc[0,1]-(sigmaln)**2))]*len(ISm.index)))    
            pisoM = pisoM/sum(pisoM)
            pisoI = pisoI/sum(pisoI)
            piso = pisoM*pisoI
            tmp.iloc[:,8] = piso/sum(piso)
            
        ### computing RT score
        pRT = [pRTNone]*len(tmp.index)
        RTranges = list(tmp['RT range'])

        for r in range(0,len(RTranges)):
            if RTranges[r]!=None:
                rtrange = RTranges[r].split(';')
                rtrange = [float(i) for i in rtrange]
                if rtm >= rtrange[0] and rtm <= rtrange[1]:
                    pRT[r]=1
                else:
                    pRT[r]=pRTout

        ### if for this id I have measured fragmentation spectrum (or spectra) and I have info in the DBMS2 for the possible hits I compute the fragementation score!
        ### add code here!
        mzid = df['ids'][k] ## getting the mzid from the dataset
        ms2inds = list(allAdds['MS2'][hits]) ## getting the MS2 DB ids 
        ### I need to do anything only if the mzid is present in the dfMS2 AND if I have and ms2inds different than None
        if mzid in list(dfMS2['id']) and (len(ms2inds)-ms2inds.count(None))>0:
            Msps = dfMS2[dfMS2['id']==mzid] ### getting the measured spectra for this id
            tmp.iloc[len(hits),9]=CSunk
            for h in range(0,len(hits)):
                #### look for the spectra in the database
                MS2id=ms2inds[h]
                CS=[]
                for s in range(0,len(Msps.index)):
                    Msp = Msps.iloc[s,1]
                    ev = Msps.iloc[s,2]
                    precType = allAdds['adduct'][hits[h]]
                    DBsp = DBMS2[(DBMS2['MonaID']==MS2id) & (DBMS2['precursorType']==precType) & (DBMS2['collision.energy']==ev)]['spectrum']
                    if len(DBsp)>0:
                        DBsp = DBsp.item()
                        CS.append(MS2compare.cosine_similarity(DBsp,Msp,mzdCS,ppmCS))
                    else:
                        CS.append(0.0)
                CS = max(CS)
                if CS==0:
                    CS=CSunk
                tmp.iloc[h,9]=CS
            

        ### computing posteriors integrating iso score and pk and RTscore and MS2 score
        p1 = tmp.iloc[:,8]
        pMS2 = tmp.iloc[:,9]
        pMS2 = [CSunk if v is None else v for v in pMS2]# fragmentation score remove none
        pMS2 = [x / sum(pMS2) for x in pMS2] # normalize fragmentation scores
        p1 = [1 if v is None else v for v in p1]# isopattern score
        pks = [1 if v is None else v for v in pks]# previous knowledge score
        post = [a*b*c*d*e for a,b,c,d,e in zip(priors,p1,pks,pRT,pMS2)]
        post = [x / sum(post) for x in post]
        tmp.iloc[:,11] = post
        tmp = tmp.sort_values(by=['post'], ascending=False)
            
    else:
        tmp=pandas.DataFrame({'id':["Unknown"],'name':["Unknown"],'formula':[None],'adduct':[None],'m/z':[None],
                       'charge':[None],'RT range':[None],'ppm':[ppm],'isotope pattern score':[None],
                       'fragmentation pattern score':[None],'prior':[1],'post':[1]})
    tmp.index=range(0,len(tmp.index)) 
    return(tmp)

def Gibbs_sampler_add(df,annotations,noits=100,burn=None,delta_add=1,all_out=False,zs=None):
    """Gibbs sampler considering only adduct connections. The function computes the posterior probabilities
    of the annotations considering the adducts connections.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
                     Output of functions MS1annotation(), MS1annotation_Parallel(), MSMSannotation() or MSMSannotation_Parallel
        noits: number of iterations if the Gibbs sampler to be run
        burn: number of iterations to be ignored when computing posterior probabilities. If None, is set to 10% of total iterations
        delta_add: parameter used when computiong the conditional priors. The parameter must be positive. The smaller the parameter
                   the more weight the adducts connections have on the posterior probabilities. Default 1.
        all_out: logical value. If true the list of assigments found in each iteration is returned by the function. Default False.
        zs: list of assigments computed in a previous run of the Gibbs sampler. Optional, default None.
    Output:
        annotations: the function modifies the annotations dictionary by adding 2 columns to each entry. One named 'post Gibbs'
                     contains the posterior probabilities computed. The other is called 'chi-square pval' containing the p-value from
                     a chi-squared test comparing the 'post' with the 'post Gibbs' probabilities.
        zs: optional, if all_out==True, the function return the full list of assignments computed. This allows restarting the sampler
            from where you are from a previous run.
    """
    start = time.time()
    print("computing posterior probabilities including adducts connections")
    print("initialising sampler ...")
    noits = int(noits) 
    
    ks = list(annotations.keys())
    rids = [] #get a vector of relation ids associated with the annotated features
    for k in ks:
        rids.append(df[df['ids']==k]['rel.ids'].item())
    ca = [] # initialise current annotation vector
    ca_id = []



    if zs is None:
        for k in ks: ### go through each mass for which I have an annotation I CAN PARALLELIZE THIS
            tmp = annotations[k]
            P=list(tmp['post'])### vector of the probabilities currently assigned to each possible annotation
            a_list = list(range(0,len(P))) ### I used this as vector of assignments
            c=random.choices(a_list, P) ### for the mass k, randomly choose an annotation based on probabilities on P
            ca.append(c) ### store the index
            ca_id.append(tmp.iloc[c,0].item()) ### store the id   
        zs = []
        zs.append(ca.copy())
        noits2=noits 
    else:
        ca = zs[len(zs)-1]
        ca_id = []
        for i in range(0,len(ca)):
            k=ks[i]
            tmp = annotations[k]
            ca_id.append(tmp.iloc[[ca[i]],0].item()) ### store the id
        noits2=noits+len(zs)
    
    indk = list(range(0,len(ks)))
    for it in tqdm(range(0,noits), desc = 'Gibbs Sampler Progress Bar'):
        ca, ca_id = gibbs_sampler_add_iter(indk,ks,rids,annotations,ca_id,ca,delta_add,it)
        zs.append(ca.copy())          
            
    zsdf = pandas.DataFrame(zs).transpose()
    
    if burn is None:
        burn = int(noits2*0.10)
    print('parsing results ...')
    positions = range(burn,noits2)                                                   
    for m in range(0,len(zsdf.index)):
        id= ks[m]
        post_gibbs = list([0.0]*annotations[id].index)
        observed = list([0.0]*annotations[id].index)
        expected = list([0.0]*annotations[id].index)
        pold= list(annotations[id]['post'])
        fs=collections.Counter(zsdf.iloc[m,positions])
        for i in fs.keys():
            post_gibbs[i] =fs[i]/len(positions)
            observed[i] = fs[i]
            expected[i] = pold[i]*len(positions)    
        annotations[id]['post Gibbs']=post_gibbs
        keep = [i for i, e in enumerate(expected) if e != 0]
        expected=[expected[i] for i in keep]
        observed = [observed[i] for i in keep]
        expected = [x+((sum(observed)-sum(expected))/len(expected)) for x in expected] ## when computing the expected frequencies there are numerical problems....
        res = stats.chisquare(f_obs=observed, f_exp=expected)
        annotations[id]['chi-square pval']= res.pvalue
        if res.pvalue < 0.001:
            annotations[id]=annotations[id].sort_values(by=['post Gibbs'], ascending=False)
    
    end = time.time()
    print('Done - ',round(end - start,1), 'seconds elapsed')
    if all_out:
        return(zs)


def gibbs_sampler_add_iter(indk,ks,rids,annotations,ca_id,ca,delta_add,it):
    random.shuffle(indk) #### I need to randomize the order I use to go through each mass in each iteration
    for i in indk: 
        k = ks[i]
        rid=rids[i]
        tmp = annotations[k] ### get the annotation table for this mass
        p = list(tmp['post']) ### get the prior probabilities from the annotation table.
        p_add = [0]*len(tmp.index) ### initialize the p_add vector
        ca_id2=[] # contains the the current annotation ids for the masses with the same the relation id
        for r in range(0,len(rids)):
            if rids[r]==rid and r!=k:
                ca_id2.append(ca_id[r]) ## 

        #padd = [ca_id2.count(x) for x in idcps]
        for cp in range(0,len(p_add)): ###populating p_add by counting add connections
            idcp = tmp.iloc[cp,0] ### annotation considered for this position of p_add
            p_add[cp] = ca_id2.count(idcp)## count the number of add connections
        p_add=[x+delta_add for x in p_add] ###computing the actual p_add 1/2
        p_add=[x/sum(p_add) for x in p_add] ###computing the actual p_add 2/2
        p0= [a * b for a, b in zip(p, p_add)] ### merging pbio and prior
        p0=[x/sum(p0) for x in p0] ## normalize merged probabilities
        a_list = list(range(0,len(p0))) ### I used this as vector of assignments
        c=random.choices(a_list, p0) ### randomly choose an annotation based on probabilities on merged prob
        ca[i] = c[0] # update current annotation
        ca_id[i] = tmp.iloc[c,0].item() # update current annotation id
    return(ca, ca_id)




def bio_single_iter_reactions(all_ids_DB,all_rs_DB,id1,id2):
    r1= all_rs_DB[all_ids_DB.index(id1)]
    r2= all_rs_DB[all_ids_DB.index(id2)]
    r1=r1.split()
    r2=r2.split()
    out = len(set.intersection(set(r1),set(r2)))
    if out >0:
        return((id1,id2))
    else:
        return(('x','x'))





def bio_single_iter_connections(all_ids_DB,all_forms_DB,conns,id1,id2):
    f1= molmass.Formula(all_forms_DB[all_ids_DB.index(id1)])
    f2= molmass.Formula(all_forms_DB[all_ids_DB.index(id2)])
    diff= None
    if util.check_ded(f1,f2):
        diff = f1.__sub__(f2).formula
    elif util.check_ded(f2,f1):
        diff = f2.__sub__(f1).formula        
    if diff in conns:
        return((id1,id2))
    else:
        return(('x','x'))

    
    
    
def Compute_Bio(DB, annotations, mode='reactions', connections = ["C3H5NO", "C6H12N4O", "C4H6N2O2", "C4H5NO3", "C3H5NOS", "C6H10N2O3S2",
               "C5H7NO3","C5H8N2O2","C2H3NO","C6H7N3O","C6H11NO","C6H11NO","C6H12N2O",
               "C5H9NOS","C9H9NO","C5H7NO","C3H5NO2","C4H7NO2","C11H10N2O","C9H9NO2",
               "C5H9NO","C4H4O2","C3H5O","C10H12N5O6P","C10H15N2O3S","C10H14N2O2S","CH2ON",
               "C21H34N7O16P3S","C21H33N7O15P3S","C10H15N3O5S","C5H7","C3H2O3","C16H30O",
               "C8H8NO5P","CH3N2O","C5H4N5","C10H11N5O3","C10H13N5O9P2","C10H12N5O6P",
               "C9H13N3O10P2","C9H12N3O7P","C4H4N3O","C10H13N5O10P2","C10H12N5O7P","C5H4N5O",
               "C10H11N5O4","C10H14N2O10P2","C10H12N2O4","C5H5N2O2","C10H13N2O7P","C9H12N2O11P2",
               "C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2","C2H2O","C2H2","CO2","CHO2","H2O","H3O6P2",
               "C2H4","CO","C2O2","H2","O","P","C2H2O","CH2","HPO3","NH2","PP","NH","SO3","N","C6H10O5",
               "C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6","C6H10O5","C18H30O15"]):
    """Compute matrix of biochemical connections. Either based on a list of possible connections in the form of a list of formulas or
    based on the reactions present in the database.
    Inputs:
        DB: pandas dataframe containing the database against which the annotation is performed. The DB must
            contain the following columns in this exact order (optional fields can contain None):
                - id: unique id of the database entry (e.g., 'C00031') - necessary
                - name: compund name (e.g., 'D-Glucose') - necessary
                - formula: chemical formula (e.g., 'C6H12O6') - necessary
                - inchi: inchi string - optional
                - smiles: smiles string - optional
                - RT: if known, retention time range (in seconds) where this compound is expected
                      to elute (e.g., '30;60') - optional
                - adducts: list of adducts that should be considered for this entry (e.g.,'M+Na;M+H;M+')
                - description: comments on the entry - optional
                - pk: previous knowledge on the likelihood of this compoud to be present in the sample
                      analyse. The value has to be between 1 (compound likely to be present in the sample)
                      and 0 (compound cannot be present in the sample).
                - MS2: id for the MS2 database entries related to this compound (optional)
                - reactions: list of reactions ids involving this compound (e.g., 'R00010 R00015 R00028')-optional, but necessary if
                  mode='reactions'.
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
                     Output of functions MS1annotation(), MS1annotation_Parallel(), MSMSannotation() or MSMSannotation_Parallel
        mode: either 'reactions' (connections are computed based on the reactions present in the database) or 'connections' (connections are
              computed based on the list of connections provided). Default 'reactions'.
        connections: list of possible connections between compouds defined as formulas. Only necessary if mode='connections'. A list of common
                     biotransformations is provided as default.
    Output:
        Bio: dataframe containing all the possible connections computed.
    """
    print("computing all possible biochemical connections")
    start = time.time()
    DB.loc[DB.reactions==None,'reactions']=''
        
    ### getting the ids of all possible hits to the database
    all_Ks = list(annotations.keys())
    all_ids = []
    for k in all_Ks:
        all_ids= all_ids+annotations[k]['id'].tolist()
    
    all_ids = list(set(all_ids))
    all_ids.remove('Unknown')
    all_ids_DB = DB['id'].to_list()
    Bio = []
    
    
    
    if mode=='connections':
        print("considering the provided connections ...")
        all_forms_DB = DB['formula'].to_list()
        conns = []
        for c in connections:
            conns.append(molmass.Formula(c).formula)
        
        
        for x,y in itertools.combinations(all_ids, 2):
            Bio = Bio+[bio_single_iter_connections(all_ids_DB,all_forms_DB,connections,x,y)]


    
    elif mode=='reactions':
        print("considering the reactions stored in the database ...")
        all_rs_DB = DB['reactions'].to_list()
        all_rs_DB= ['' if v is None else v for v in all_rs_DB]
        for x,y in itertools.combinations(all_ids, 2):
            Bio = Bio+[bio_single_iter_reactions(all_ids_DB,all_rs_DB,x,y)]


    
    else:
        print('ERROR: mode can only be connections or reactions')
        return

    Bio = [i for i in Bio if i != ('x','x')]
    Bio = pandas.DataFrame(Bio, index=None)
    end = time.time()   
    print(round(end - start,1), 'seconds elapsed')
    return(Bio)


def Compute_Bio_Parallel(DB, annotations, mode='reactions', connections = ["C3H5NO", "C6H12N4O", "C4H6N2O2", "C4H5NO3", "C3H5NOS", "C6H10N2O3S2",
               "C5H7NO3","C5H8N2O2","C2H3NO","C6H7N3O","C6H11NO","C6H11NO","C6H12N2O",
               "C5H9NOS","C9H9NO","C5H7NO","C3H5NO2","C4H7NO2","C11H10N2O","C9H9NO2",
               "C5H9NO","C4H4O2","C3H5O","C10H12N5O6P","C10H15N2O3S","C10H14N2O2S","CH2ON",
               "C21H34N7O16P3S","C21H33N7O15P3S","C10H15N3O5S","C5H7","C3H2O3","C16H30O",
               "C8H8NO5P","CH3N2O","C5H4N5","C10H11N5O3","C10H13N5O9P2","C10H12N5O6P",
               "C9H13N3O10P2","C9H12N3O7P","C4H4N3O","C10H13N5O10P2","C10H12N5O7P","C5H4N5O",
               "C10H11N5O4","C10H14N2O10P2","C10H12N2O4","C5H5N2O2","C10H13N2O7P","C9H12N2O11P2",
               "C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2","C2H2O","C2H2","CO2","CHO2","H2O","H3O6P2",
               "C2H4","CO","C2O2","H2","O","P","C2H2O","CH2","HPO3","NH2","PP","NH","SO3","N","C6H10O5",
               "C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6","C6H10O5","C18H30O15"], ncores=1):
    """Compute possible biochemical connetions given a database and a list of
    possible connectios in the form of a list of formulas
    """
    print("computing all possible biochemical connections")
    start = time.time()
    DB.loc[DB.reactions==None,'reactions']=''

    ### getting the ids of all possible hits to the database
    all_Ks = list(annotations.keys())
    all_ids = []
    for k in all_Ks:
        all_ids= all_ids+annotations[k]['id'].tolist()
    
    all_ids = list(set(all_ids))
    all_ids.remove('Unknown')
    all_ids_DB = DB['id'].to_list()
    
    
    
    if mode=='connections':
        print("considering the provided connections ...")
        all_forms_DB = DB['formula'].to_list()
        conns = []
        for c in connections:
            conns.append(molmass.Formula(c).formula)
        
        pool_obj = multiprocessing.Pool(ncores)
        Bio = pool_obj.starmap(partial(bio_single_iter_connections,all_ids_DB,all_forms_DB,connections),itertools.combinations(all_ids, 2))
        pool_obj.terminate()

    
    elif mode=='reactions':
        print("considering the reactions stored in the database ...")
        all_rs_DB = DB['reactions'].to_list()
        all_rs_DB= ['' if v is None else v for v in all_rs_DB]
        pool_obj = multiprocessing.Pool(ncores)
        Bio = pool_obj.starmap(partial(bio_single_iter_reactions,all_ids_DB,all_rs_DB),itertools.combinations(all_ids, 2))
        pool_obj.terminate()

    
    else:
        print('ERROR: mode can only be connections or reactions')
        return

    
    Bio = [i for i in Bio if i != ('x','x')]
    Bio = pandas.DataFrame(Bio, index=None)
    end = time.time()   
    print(round(end - start,1), 'seconds elapsed')
    return(Bio)



def Gibbs_sampler_bio(df,annotations,Bio,noits=100,burn=None,delta_bio=1,all_out=False,zs=None):
    """Gibbs sampler considering only possible biochemical connections. The function computes the
    posterior probabilities of the annotations considering the possible biochemical connections
    reported in Bio.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
                     Output of functions MS1annotation(), MS1annotation_Parallel(), MSMSannotation() or MSMSannotation_Parallel
        Bio: dataframe (2 columns), reporting all the possible connections between compounds. It uses the unique ids from the database.
             It could be the output of Compute_Bio() or Compute_Bio_Parallel().
        noits: number of iterations if the Gibbs sampler to be run
        burn: number of iterations to be ignored when computing posterior probabilities. If None, is set to 10% of total iterations
        delta_bio: parameter used when computiong the conditional priors. The parameter must be positive. The smaller the parameter
                   the more weight the adducts connections have on the posterior probabilities. Default 1.
        all_out: logical value. If true the list of assigments found in each iteration is returned by the function. Default False.
        zs: list of assigments computed in a previous run of the Gibbs sampler. Optional, default None.
    Output:
        annotations: the function modifies the annotations dictionary by adding 2 columns to each entry. One named 'post Gibbs'
                     contains the posterior probabilities computed. The other is called 'chi-square pval' containing the p-value from
                     a chi-squared test comparing the 'post' with the 'post Gibbs' probabilities.
        zs: optional, if all_out==True, the function return the full list of assignments computed. This allows restarting the sampler
            from where you are from a previous run.
    """
    start = time.time()
    print("computing posterior probabilities including biochemical connections")
    print("initialising sampler ...")
    noits = int(noits) 
    
    ks = list(annotations.keys())
    rids = [] #get a vector of relation ids associated with the annotated features
    for k in ks:
        rids.append(df[df['ids']==k]['rel.ids'].item())
    ca = [] # initialise current annotation vector
    ca_id = []
    Bio = list(Bio.itertuples(index=False, name=None))


    if zs is None:
        for k in ks: ### go through each mass for which I have an annotation I CAN PARALLELIZE THIS
            tmp = annotations[k]
            P=list(tmp['post']) ### vector of the probabilities currently assigned to each possible annotation
            a_list = list(range(0,len(P))) ### I used this as vector of assignments
            c=random.choices(a_list, P) ### for the mass k, randomly choose an annotation based on probabilities on P
            ca.append(c) ### store the index
            ca_id.append(tmp.iloc[c,0].item()) ### store the id   
        zs = []
        zs.append(ca.copy())
        noits2=noits 
    else:
        ca = zs[len(zs)-1]
        ca_id = []
        for i in range(0,len(ca)):
            k=ks[i]
            tmp = annotations[k]
            ca_id.append(tmp.iloc[[ca[i]],0].item()) ### store the id
        noits2=noits+len(zs)
    
    indk = list(range(0,len(ks)))
    for it in tqdm(range(0,noits), desc = 'Gibbs Sampler Progress Bar'):
        ca, ca_id = gibbs_sampler_bio_iter(indk,ks,annotations,Bio,ca_id,ca,delta_bio,it)
        zs.append(ca.copy())          
            
    zsdf = pandas.DataFrame(zs).transpose()
    
    if burn is None:
        burn = int(noits2*0.10)

    print('parsing results ...')
    positions = range(burn,noits2)                                                     
    for m in range(0,len(zsdf.index)):
        id= ks[m]
        post_gibbs = list([0.0]*annotations[id].index)
        observed = list([0.0]*annotations[id].index)
        expected = list([0.0]*annotations[id].index)
        pold= list(annotations[id]['post'])
        fs=collections.Counter(zsdf.iloc[m,positions])
        for i in fs.keys():
            post_gibbs[i] =fs[i]/len(positions)
            observed[i] = fs[i]
            expected[i] = pold[i]*len(positions)    
        annotations[id]['post Gibbs']=post_gibbs
        keep = [i for i, e in enumerate(expected) if e != 0]
        expected=[expected[i] for i in keep]
        observed = [observed[i] for i in keep]
        expected = [x+((sum(observed)-sum(expected))/len(expected)) for x in expected] ## when computing the expected frequencies there are numerical problems....
        res = stats.chisquare(f_obs=observed, f_exp=expected)
        annotations[id]['chi-square pval']= res.pvalue
        if res.pvalue < 0.001:
            annotations[id]=annotations[id].sort_values(by=['post Gibbs'], ascending=False)
    
    end = time.time()
    print('Done - ',round(end - start,1), 'seconds elapsed')
    if all_out:
        return(zs)



def gibbs_sampler_bio_iter(indk,ks,annotations,Bio,ca_id,ca,delta_bio,it):
    random.shuffle(indk)
    for i in indk:
        k=ks[i]
        tmp = annotations[k]
        p = list(tmp['post']) ### get the prior probabilities from the annotation table.
        p_bio = [0]*len(tmp.index) ### initialize the p_bio vector
        ca_id2=ca_id.copy() # contains the the current annotation
        del ca_id2[i] # without the annotation for the mass considered 

        for cp in range(0,len(p_bio)): ###populating p_bio by counting bio connections
            idcp = tmp.iloc[cp,0] ### annotation considered for this position of p_bio
            idcp_list = [idcp]*len(ca_id2)
            A = list(zip(idcp_list,ca_id2))+list(zip(ca_id2,idcp_list))
            p_bio[cp] = len(set.intersection(set(Bio),set(A)))
        p_bio=[x+delta_bio for x in p_bio] ###computing the actual p_bio 1/2
        p_bio=[x/sum(p_bio) for x in p_bio] ###computing the actual p_bio 2/2
        p0= [a * b for a, b in zip(p, p_bio)] ### merging pbio and prior
        p0=[x/sum(p0) for x in p0] ## normalize merged probabilities
        a_list = list(range(0,len(p0))) ### I used this as vector of assignments
        c=random.choices(a_list, p0) ### randomly choose an annotation based on probabilities on merged prob
        ca[i] = c[0] # update current annotation
        ca_id[i] = tmp.iloc[c,0].item() # update current annotation id
    return(ca, ca_id)



def Gibbs_sampler_bio_add(df,annotations,Bio,noits=100,burn=None,delta_bio=1,delta_add=1,all_out=False,zs=None):
    """Gibbs sampler considering both biochemical and adducts connections. The function computes the
    posterior probabilities of the annotations considering the possible biochemical connections
    reported in Bio and the possible adducts connection.
    Inputs:
        df: pandas dataframe containg the MS1 data. It should be the output of the function ipa.map_isotope_patterns()
        annotations: a dictonary containg all the possible annotations for the measured features. The keys of the dictionay are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
                     Output of functions MS1annotation(), MS1annotation_Parallel(), MSMSannotation() or MSMSannotation_Parallel
        Bio: dataframe (2 columns), reporting all the possible connections between compounds. It uses the unique ids from the database.
             It could be the output of Compute_Bio() or Compute_Bio_Parallel().
        noits: number of iterations if the Gibbs sampler to be run
        burn: number of iterations to be ignored when computing posterior probabilities. If None, is set to 10% of total iterations
        delta_bio: parameter used when computiong the conditional priors. The parameter must be positive. The smaller the parameter
                   the more weight the adducts connections have on the posterior probabilities. Default 1.
        delta_add: parameter used when computiong the conditional priors. The parameter must be positive. The smaller the parameter
                   the more weight the biochemical connections have on the posterior probabilities. Default 1.
        all_out: logical value. If true the list of assigments found in each iteration is returned by the function. Default False.
        zs: list of assigments computed in a previous run of the Gibbs sampler. Optional, default None.
    Output:
        annotations: the function modifies the annotations dictionary by adding 2 columns to each entry. One named 'post Gibbs'
                     contains the posterior probabilities computed. The other is called 'chi-square pval' containing the p-value from
                     a chi-squared test comparing the 'post' with the 'post Gibbs' probabilities.
        zs: optional, if all_out==True, the function return the full list of assignments computed. This allows restarting the sampler
            from where you are from a previous run.
    """
    start = time.time()
    print("computing posterior probabilities including biochemical and adducts connections")
    print("initialising sampler ...")
    noits = int(noits) 
    ks = list(annotations.keys())
    rids = [] #get a vector of relation ids associated with the annotated features
    for k in ks:
        rids.append(df[df['ids']==k]['rel.ids'].item())
    ca = [] # initialise current annotation vector
    ca_id = []
    Bio = list(Bio.itertuples(index=False, name=None))


    if zs is None:
        for k in ks: ### go through each mass for which I have an annotation I CAN PARALLELIZE THIS
            tmp = annotations[k]
            P=list(tmp['post']) ### vector of the probabilities currently assigned to each possible annotation
            a_list = list(range(0,len(P))) ### I used this as vector of assignments
            c=random.choices(a_list, P) ### for the mass k, randomly choose an annotation based on probabilities on P
            ca.append(c) ### store the index
            ca_id.append(tmp.iloc[c,0].item()) ### store the id   
        zs = []
        zs.append(ca.copy())
        noits2=noits 
    else:
        ca = zs[len(zs)-1]
        ca_id = []
        for i in range(0,len(ca)):
            k=ks[i]
            tmp = annotations[k]
            ca_id.append(tmp.iloc[[ca[i]],0].item()) ### store the id
        noits2=noits+len(zs)
    
    indk = list(range(0,len(ks)))
    for it in tqdm(range(0,noits), desc = 'Gibbs Sampler Progress Bar'):
        ca, ca_id = gibbs_sampler_bio_add_iter(indk,ks,rids,annotations,Bio,ca_id,ca,delta_bio,delta_add,it)
        zs.append(ca.copy())          
            
    zsdf = pandas.DataFrame(zs).transpose()
    
    if burn is None:
        burn = int(noits2*0.10)
    print('parsing results ...')
    positions = range(burn,noits2)                                                     
    for m in range(0,len(zsdf.index)):
        id= ks[m]
        post_gibbs = list([0.0]*annotations[id].index)
        observed = list([0.0]*annotations[id].index)
        expected = list([0.0]*annotations[id].index)
        pold= list(annotations[id]['post'])
        fs=collections.Counter(zsdf.iloc[m,positions])
        for i in fs.keys():
            post_gibbs[i] =fs[i]/len(positions)
            observed[i] = fs[i]
            expected[i] = pold[i]*len(positions)    
        annotations[id]['post Gibbs']=post_gibbs
        keep = [i for i, e in enumerate(expected) if e != 0]
        expected=[expected[i] for i in keep]
        observed = [observed[i] for i in keep]
        expected = [x+((sum(observed)-sum(expected))/len(expected)) for x in expected] ## when computing the expected frequencies there are numerical problems....
        res = stats.chisquare(f_obs=observed, f_exp=expected)
        annotations[id]['chi-square pval']= res.pvalue
        if res.pvalue < 0.001:
            annotations[id]=annotations[id].sort_values(by=['post Gibbs'], ascending=False)

    
    end = time.time()
    print('Done - ',round(end - start,1), 'seconds elapsed')
    if all_out:
        return(zs)

def gibbs_sampler_bio_add_iter(indk,ks,rids,annotations,Bio,ca_id,ca,delta_bio,delta_add,it):
    random.shuffle(indk) #### I need to randomize the order I use to go through each mass in each iteration
    for i in indk: 
        k = ks[i]
        rid=rids[i]
        tmp = annotations[k] ### get the annotation table for this mass
        p = list(tmp['post']) ### get the prior probabilities from the annotation table.
        p_bio = [0]*len(tmp.index) ### initialize the p_add vector
        p_add = [0]*len(tmp.index) ### initialize the p_add vector
        ca_id2=[] # contains the the current annotation ids for the masses with the same the relation id
        for r in range(0,len(rids)):
            if rids[r]==rid and r!=k:
                ca_id2.append(ca_id[r]) ## 

        #padd = [ca_id2.count(x) for x in idcps]
        for cp in range(0,len(p_add)): ###populating p_add by counting add connections
            idcp = tmp.iloc[cp,0] ### annotation considered for this position of p_add
            p_add[cp] = ca_id2.count(idcp)## count the number of add connections
            idcp_list = [idcp]*len(ca_id2)
            A = list(zip(idcp_list,ca_id2))+list(zip(ca_id2,idcp_list))
            p_bio[cp] = len(set.intersection(set(Bio),set(A)))
        p_add=[x+delta_add for x in p_add] ###computing the actual p_add 1/2
        p_add=[x/sum(p_add) for x in p_add] ###computing the actual p_add 2/2
        p_bio=[x+delta_bio for x in p_bio] ###computing the actual p_bio 1/2
        p_bio=[x/sum(p_bio) for x in p_bio] ###computing the actual p_bio 2/2
        p0= [a * b * c for a, b, c in zip(p, p_add, p_bio)] ### merging pbio and prior
        p0=[x/sum(p0) for x in p0] ## normalize merged probabilities
        a_list = list(range(0,len(p0))) ### I used this as vector of assignments
        c=random.choices(a_list, p0) ### randomly choose an annotation based on probabilities on merged prob
        ca[i] = c[0] # update current annotation
        ca_id[i] = tmp.iloc[c,0].item() # update current annotation id
    return(ca, ca_id)


def simpleIPA(df,ionisation,DB,adducts_file,ppm,dfMS2=None,DBMS2=None,noits=100,burn=None,delta_add=None,delta_bio=None,Bio=None,
mode='reactions',CSunk=0.5,isodiff=1,ppmiso=100,ncores=1,me=5.48579909065e-04,ratiosd=0.9,ppmunk=None,ratiounk=None,ppmthr=None,
pRTNone=None,pRTout=None,mzdCS=0, ppmCS=10,connections = ["C3H5NO", "C6H12N4O", "C4H6N2O2", "C4H5NO3", "C3H5NOS", "C6H10N2O3S2",
               "C5H7NO3","C5H8N2O2","C2H3NO","C6H7N3O","C6H11NO","C6H11NO","C6H12N2O",
               "C5H9NOS","C9H9NO","C5H7NO","C3H5NO2","C4H7NO2","C11H10N2O","C9H9NO2",
               "C5H9NO","C4H4O2","C3H5O","C10H12N5O6P","C10H15N2O3S","C10H14N2O2S","CH2ON",
               "C21H34N7O16P3S","C21H33N7O15P3S","C10H15N3O5S","C5H7","C3H2O3","C16H30O",
               "C8H8NO5P","CH3N2O","C5H4N5","C10H11N5O3","C10H13N5O9P2","C10H12N5O6P",
               "C9H13N3O10P2","C9H12N3O7P","C4H4N3O","C10H13N5O10P2","C10H12N5O7P","C5H4N5O",
               "C10H11N5O4","C10H14N2O10P2","C10H12N2O4","C5H5N2O2","C10H13N2O7P","C9H12N2O11P2",
               "C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2","C2H2O","C2H2","CO2","CHO2","H2O","H3O6P2",
               "C2H4","CO","C2O2","H2","O","P","C2H2O","CH2","HPO3","NH2","PP","NH","SO3","N","C6H10O5",
               "C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6","C6H10O5","C18H30O15"]):
    # mapping isotopes
    map_isotope_patterns(df,isoDiff=isodiff, ppm=ppmiso,ionisation=ionisation)
    # computing all adducts
    if ncores>1:
        allAdds = compute_all_adducts_Parallel(adducts_file= adducts_file, DB=DB,ncores=ncores)
    else:
        allAdds = compute_all_adducts(adducts_file= adducts_file, DB=DB)

    # computing priors
    if (dfMS2 is None) or (DBMS2 is None):
        if ncores>1:
            annotations = MS1annotation_Parallel(df=df,allAdds=allAdds,ppm=ppm,me=me,ratiosd=ratiosd,ppmunk=ppmunk,
            ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout,ncores=ncores)
        else:
            annotations = MS1annotation(df=df,allAdds=allAdds,ppm=ppm,me=me,ratiosd=ratiosd,ppmunk=ppmunk,
            ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout)
    else:
        if ncores>1:
            annotations = MSMSannotation_Parallel(df=df,dfMS2=dfMS2,allAdds=allAdds,DBMS2=DBMS2,ppm=ppm,me=me,ratiosd=ratiosd,
            ppmunk=ppmunk,ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout,mzdCS=mzdCS,ppmCS=ppmCS,
            CSunk=CSunk,ncores=ncores)
        else:
            annotations = MSMSannotation(df=df,dfMS2=dfMS2,allAdds=allAdds,DBMS2=DBMS2,ppm=ppm,me=me,ratiosd=ratiosd,
            ppmunk=ppmunk,ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout,mzdCS=mzdCS,ppmCS=ppmCS,
            CSunk=CSunk)

    # computing Bio matrix (if necessary)
    if (Bio is None) and (delta_bio is not None):
        if ncores>1:
            Bio=Compute_Bio_Parallel(DB=DB,annotations=annotations,mode=mode,connections=connections,ncores=ncores)
        else:
            Bio=Compute_Bio(DB=DB,annotations=annotations,mode=mode,connections=connections)

    # Gibbs sampler (if needed). Which one based on the inputs
    if (Bio is not None) and (delta_bio is not None) and (delta_add is not None):
        Gibbs_sampler_bio_add(df=df,annotations=annotations,Bio=Bio,noits=noits,burn=burn,delta_bio=delta_bio,delta_add=delta_add)
    elif (Bio is not None) and (delta_bio is not None) and (delta_add is None):
        Gibbs_sampler_bio(df=df,annotations=annotations,Bio=Bio,noits=noits,burn=burn,delta_bio=delta_bio)
    elif (Bio is None) and (delta_bio is None) and (delta_add is not None):
        Gibbs_sampler_add(df=df,annotations=annotations,noits=noits,burn=burn,delta_add=delta_add)
    
    return(annotations)

        

