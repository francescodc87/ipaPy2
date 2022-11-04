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
from tqdm import tqdm
from ipaPy2 import util
from ipaPy2 import iterations

__author__ = "Francesco Del Carratore"
__maintainer__ = "Francesco Del Carratore"
__email__ = "francesco.delcarratore@gmail.com"

def clusterFeatures(df,Cthr=0.8,RTwin=1,Intmode='max'):
    """
    Clustering MS1 features based on correlation across samples.
    
    Parameters
    ----------
    df: pandas dataframe with the following columns:
        -ids: a unique id for each feature
        -mzs: mass-to-charge ratios, usually the average across different
              samples.
        -RTs: retention times in seconds, usually the average across different
              samples.
        -Intensities: for each sample, a column reporting the detected
                      intensities in each sample. 
    Cthr: Default value 0.8. Minimum correlation allowed in each cluster
    RTwin: Default value 1. Maximum difference in RT time between features in
           the same cluster
    Intmode: Defines how the representative intensity of each feature is
             computed. If 'max' (default) the maximum across samples is used.
             If 'ave' the average across samples is computed
    Returns
    -------
    df: pandas dataframe in correct format to be used as an input of the
    map_isotope_patterns() function
    """
    print("Clustering features ....")
    start = time.time()
    df=df.replace('None',None)
    ids = list(df.iloc[:,0])
    mzs=list(df.iloc[:,1])
    RTs=list(df.iloc[:,2])
    CorrInts = df.iloc[:,3:len(df.index)].transpose().corr()
    fids=[]
    rids=[]
    fmz=[]
    fRTs =[]
    fInt=[]
    
    if Intmode=='max':
        Int=df.iloc[:,3:len(df.index)].max(axis=1)
    elif Intmode=='ave':
        Int=df.iloc[:,3:len(df.index)].mean(axis=1)
    else:
        raise ValueError("Intmode not allowed")
    flag = True
    rid = 0
    while(flag):
        RTdiff = [abs(v-RTs[0]) for v in RTs]
        ind =[v for v in range(0,len(RTdiff)) if (CorrInts.iloc[v,0]>=Cthr and RTdiff[v]<=RTwin)]
        indids=[ids[v] for v in ind]
        fids.extend(indids)
        rids.extend([rid]*len(indids))
        fmz.extend([mzs[v] for v in ind])
        fRTs.extend([RTs[v] for v in ind])
        fInt.extend([Int[v] for v in ind])
        rid=rid+1
        indkeep = [v for v in range(0,len(ids)) if ids[v] not in indids]
        if len(indkeep)>0:
            ids=[ids[v] for v in indkeep]
            mzs=[mzs[v] for v in indkeep]
            RTs=[RTs[v] for v in indkeep]
            Int=[Int[v] for v in indkeep]
            CorrInts=CorrInts.iloc[indkeep,indkeep]
        else:
            flag=False
        
    df = pandas.DataFrame(list(zip(fids,rids,fmz,fRTs,fInt)),
                              columns=['ids','rel.ids','mzs','RTs','Int'])
    end = time.time()
    print(round(end - start,1), 'seconds elapsed')
    return(df)



def map_isotope_patterns(df,isoDiff=1, ppm=100, ionisation=1):
    """
    mapping isotope patterns in MS1 data.
    
    Parameters
    ----------
    df : pandas dataframe (necessary)
         A dataframe containing the MS1 data including the following columns:
            -ids: an unique id for each feature
            -rel.ids:   relation ids. In a previous step of the data processing
                        pipeline, features are clustered based on peak shape
                        similarity/retention time. Features in the same
                        cluster are likely to come from the same metabolite.
                        All isotope patterns must be in the same rel.id
                        cluster.
            -mzs: mass-to-charge ratios, usually the average across
                  different samples.
            -RTs: retention times in seconds, usually the average across
                  different samples.
            -Ints: representative (e.g., maximum or average) intensity detected
                   for each feature across samples (either peak area or peak
                   intensity)
    isoDiff : Default value 1. Difference between isotopes of charge 1, does
              not need to be exact
    ppm:   Default value 100. Maximum ppm value allowed between 2 isotopes.
            It is very high on purpose
    ionisation: Default value 1. positive = 1, negative = -1
    
    Returns
    -------
    df: the main input is modified by adding and populating the following
        columns
        - relationship: the possible values are:
                        * bp: basepeak, most intense peak within each rel id
                        * bp|isotope: isotope of the basepeak
                        * potential bp: most intense peak within each isotope
                                        pattern (excluding the basepeak)
                        * potential bp|isotope: isotope of one potential bp
        - isotope pattern: feature used to cluster the different isotope
                            patterns within the same relation id
        - charge: predicted charge based on the isotope pattern (1,2,3,4,5 or
                  -1,-2,-3,-4,-5 are the only values allowed)
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
                ppm4 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/4)))/(dfg.iloc[k,2]+(isoDiff/4)))*(10**6))
                ppm5 = abs(((dfg.iloc[0:(len(dfg.index)),2]-(dfg.iloc[k,2]+(isoDiff/5)))/(dfg.iloc[k,2]+(isoDiff/5)))*(10**6))
                indiso1 = util.which(ppm1 <= ppm)
                if k in indiso1: indiso1.remove(k)
                indiso2 = util.which(ppm2 <= ppm)
                if k in indiso2: indiso2.remove(k)
                indiso3 = util.which(ppm3 <= ppm)
                if k in indiso3: indiso3.remove(k)
                indiso4 = util.which(ppm4 <= ppm)
                if k in indiso4: indiso4.remove(k)
                indiso5 = util.which(ppm5 <= ppm)
                if k in indiso5: indiso5.remove(k)
                if len(indiso5) >0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso5,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso5,6] = c
                    dfg.iloc[k,7] = 5*ionisation
                    dfg.iloc[indiso5,7] = 5*ionisation
                    f2=True
                elif len(indiso4) > 0:
                    dfg.iloc[k,5] = "isotope"
                    dfg.iloc[indiso4,5] = "isotope"
                    dfg.iloc[k,6] = c
                    dfg.iloc[indiso4,6] = c
                    dfg.iloc[k,7] = 4*ionisation
                    dfg.iloc[indiso4,7] = 4*ionisation
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


def compute_all_adducts(adductsAll, DB, ionisation=1, ncores=1):
    """
    compute all adducts table based on the information present in the database
    
    Parameters
    ----------
    adductsAll : pandas dataframe (necessary)
                 Dataframe containing information on all possible
                 adducts. The file must be in the same format as the example
                 provided in the DB/adducts.csv
    DB : pandas dataframe (necessary)
         Dataframe containing the database against which the annotation is
         performed. The DB must contain the following columns in this exact
         order (optional fields can contain None):
             - id: unique id of the database entry (e.g., 'C00031') - necessary
             - name: compound name (e.g., 'D-Glucose') - necessary
             - formula: chemical formula (e.g., 'C6H12O6') - necessary
             - inchi: inchi string - optional
             - smiles: smiles string - optional
             - RT: if known, retention time range (in seconds) where this
                   compound is expected to elute (e.g., '30;60') - optional
             - adductsPos: list of adducts that should be considered in
                           positive mode for this entry (e.g.,'M+Na;M+H;M+')
             - adductsNeg: list of adducts that should be considered in
                           negative mode for this entry (e.g.,'M-H;M-2H')
             - description: comments on the entry - optional
             - pk: previous knowledge on the likelihood of this compound to be
                   present in the sample analysed. The value has to be between
                   1 (compound likely to be present in the sample) and 0
                   (compound cannot be present in the sample).
             - MS2: id for the MS2 database entries related to this compound
                    (optional)
             - reactions: list of reactions ids involving this compound
                          (e.g., 'R00010 R00015 R00028')-optional 
    ionisation : Default value 1. positive = 1, negative = -1
    ncores : default value 1. Number of cores used
    
    Returns
    -------
    allAdds: pandas dataframe containing the information on all the possible
    adducts given the database.
    """
    if ncores==1:
        print("computing all adducts ....")
        start = time.time()
        DB = DB.replace(numpy.nan,None)
        data=[]
        for db in range(0,len(DB.index)):
            data.append(iterations.all_adducts_iter(DB,adductsAll,ionisation,db))
        allAdds =pandas.concat(data,ignore_index=True)
        allAdds.columns=['id','name','adduct','formula','charge','m/z','RT','pk','MS2']
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    elif ncores>1:
        print("computing all adducts - Parallelized ....")
        start = time.time()
        DB = DB.replace(numpy.nan,None)
        pool_obj = multiprocessing.Pool(ncores)
        data = pool_obj.map(partial(iterations.all_adducts_iter,DB,adductsAll,ionisation),range(0,len(DB.index)))
        pool_obj.terminate()
        allAdds =pandas.concat(data,ignore_index=True)
        allAdds.columns=['id','name','adduct','formula','charge','m/z','RT','pk','MS2']
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    else:
        raise ValueError("ncores must be >=1")
    return(allAdds)


 

def MS1annotation(df,allAdds,ppm,me = 5.48579909065e-04,ratiosd=0.9,
                  ppmunk=None,ratiounk=None,ppmthr=None, pRTNone=None,
                  pRTout=None,ncores=1):
    """
    Annotation of the dataset base on the MS1 information. Prior probabilities
    are based on mass only, while post probabilities are based on mass, RT,
    previous knowledge and isotope patterns.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
    allAdds: pandas dataframe containing the information on all the possible
            adducts given the database. It should be the output of either
            ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
    ppm: accuracy of the MS instrument used
    me: accurate mass of the electron. Default 5.48579909065e-04
    ratiosd: default 0.9. It represents the acceptable ratio between predicted
             intensity and observed intensity of isotopes. It is used to compute
             the shape parameters of the lognormal distribution used to
             calculate the isotope pattern scores as sqrt(1/ratiosd)
    ppmunk: ppm associated to the 'unknown' annotation. If not provided equal
            to ppm.
    ratiounk: isotope ratio associated to the 'unknown' annotation. If not
              provided equal to 0.5
    ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to
            2*ppm
    pRTNone: Multiplicative factor for the RT if no RTrange present in the
             database. If not provided equal to 0.8
    pRTout: Multiplicative factor for the RT if measured RT is outside the
            RTrange present in the database. If not provided equal to 0.4
    ncores: default value 1. Number of cores used
    
    Returns
    -------
    annotations: a dictionary containing all the possible annotations for the
                measured features. The keys of the dictionary are the unique
                ids for the features present in df. For each feature, the
                annotations are summarized in a pandas dataframe.
    """
    df=df.replace('None',None)
    if ncores==1:
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
            data.append(iterations.MS1_ann_iter(df,allAdds,ppm,me,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,sigmaln,k))
        keys = list(df.iloc[ind,0])
        annotations = dict(zip(keys, data))
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    elif ncores>1:
        print("annotating based on MS1 information - Parallelized ...")
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
        data = pool_obj.map(partial(iterations.MS1_ann_iter,df,allAdds,ppm,me,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,sigmaln),ind)
        pool_obj.terminate()
        keys = list(df.iloc[ind,0])
        annotations = dict(zip(keys, data))
        ##convert data into dictionary! annotations[df.iloc[k,0]]=tmp
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    else:
        raise ValueError("ncores must be >=1")
    return(annotations)        



def MSMSannotation(df,dfMS2,allAdds,DBMS2,ppm,me = 5.48579909065e-04,
                   ratiosd=0.9,ppmunk=None, ratiounk=None,ppmthr=None,
                   pRTNone=None, pRTout=None,mzdCS=0, ppmCS=10, CSunk=0.7,
                   evfilt=False,ncores=1):
    """
    Annotation of the dataset base on the MS1 and MS2 information. Prior
    probabilities are based on mass only, while post probabilities are based
    on mass, RT, previous knowledge and isotope patterns.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
    dfMS2: pandas dataframe containing the MS2 data. It must contain 3 columns
        -id: an unique id for each feature for which the MS2 spectrum was
             acquired (same as in df)
        -spectrum: string containing the spectrum information in the following
                   format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'
        -ev: collision energy used to acquire the fragmentation spectrum
    allAdds: pandas dataframe containing the information on all the possible
            adducts given the database. It should be the output of either
            ipa.compute_all_adducts() or ipa.compute_all_adducts_Parallel()
    DBMS2: pandas dataframe containing the database containing the MS2
           information
    ppm: accuracy of the MS instrument used
    me: accurate mass of the electron. Default 5.48579909065e-04
    ratiosd: default 0.9. It represents the acceptable ratio between predicted
            intensity and observed intensity of isotopes. it is used to compute
            the shape parameters of the lognormal distribution used to
            calculate the isotope pattern scores as sqrt(1/ratiosd)
    ppmunk: ppm associated to the 'unknown' annotation. If not provided equal
            to ppm.
    ratiounk: isotope ratio associated to the 'unknown' annotation. If not
              provided equal to 0.5
    ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to
            2*ppm
    pRTNone: Multiplicative factor for the RT if no RTrange present in the
            database. If not provided equal to 0.8
    pRTout: Multiplicative factor for the RT if measured RT is outside the
            RTrange present in the database. If not provided equal to 0.4
    mzdCS: maximum mz difference allowed when computing cosine similarity
           scores. If one wants to use this parameter instead of ppmCS, this
           must be set to 0. Default 0.
    ppmCS: maximum ppm allowed when computing cosine similarity scores.
           If one wants to use this parameter instead of mzdCS, this must be
           set to 0. Default 10.
    CSunk: cosine similarity score associated with the 'unknown' annotation.
            Default 0.7
    evfilt: Default value False. If true, only spectrum acquired with the same
            collision energy are considered.
    ncores: default value 1. Number of cores used
    
    Returns
    -------
    annotations: a dictionary containing all the possible annotations for the
                 measured features. The keys of the dictionary are the unique
                 ids for the features present in df. For each feature, the
                 annotations are summarized in a pandas dataframe.
    """
    df=df.replace('None',None)
    dfMS2=dfMS2.replace('None',None)
    if ncores==1:
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
            if evfilt:
                data.append(iterations.MSMS_ann_iter1(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k))
            else:
                data.append(iterations.MSMS_ann_iter2(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k))
        keys = list(df.iloc[ind,0])
        annotations = dict(zip(keys, data))
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    elif ncores>1:
        print("annotating based on MS1 and MS2 information - Parallelized...")
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
        if evfilt:
            data = pool_obj.map(partial(iterations.MSMS_ann_iter1,df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln),ind)
        else:
            data = pool_obj.map(partial(iterations.MSMS_ann_iter2,df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln),ind)
        pool_obj.terminate()
        keys = list(df.iloc[ind,0])
        annotations = dict(zip(keys, data))
        end = time.time()
        print(round(end - start,1), 'seconds elapsed')
    else:
        raise ValueError("ncores must be >=1")

    return(annotations)       



def Gibbs_sampler_add(df,annotations,noits=100,burn=None,delta_add=1,
                      all_out=False,zs=None):
    """
    Gibbs sampler considering only adduct connections. The function computes
    the posterior probabilities of the annotations considering the adducts
    connections.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
    annotations: a dictionary containing all the possible annotations for the
                measured features. The keys of the dictionary are the unique
                ids for the features present in df. For each feature, the
                annotations are summarized in a pandas dataframe. Output of
                functions MS1annotation(), MS1annotation_Parallel(),
                MSMSannotation() or MSMSannotation_Parallel
    noits: number of iterations if the Gibbs sampler to be run
    burn: number of iterations to be ignored when computing posterior
          probabilities. If None, is set to 10% of total iterations
    delta_add: parameter used when computing the conditional priors. The
               parameter must be positive. The smaller the parameter the more
               weight the adducts connections have on the posterior
               probabilities. Default 1.
    all_out: logical value. If true the list of assignments found in each
             iteration is returned by the function. Default False.
    zs: list of assignments computed in a previous run of the Gibbs sampler. 
        Optional, default None.
    
    Returns
    -------
    annotations: the function modifies the annotations dictionary by adding 2
                 columns to each entry. One named 'post Gibbs' contains the
                 posterior probabilities computed. The other is called
                 'chi-square pval' containing the p-value from a chi-squared
                 test comparing the 'post' with the 'post Gibbs' probabilities.
    zs: optional, if all_out==True, the function return the full list of
        assignments computed. This allows restarting the sampler from where
        you are from a previous run.
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
        ca, ca_id = iterations.gibbs_sampler_add_iter(indk,ks,rids,annotations,ca_id,ca,delta_add,it)
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

    
def Compute_Bio(DB, annotations=None, mode='reactions', connections = ["C3H5NO",
               "C6H12N4O", "C4H6N2O2", "C4H5NO3", "C3H5NOS", "C6H10N2O3S2",
               "C5H7NO3","C5H8N2O2","C2H3NO","C6H7N3O","C6H11NO","C6H11NO",
               "C6H12N2O","C5H9NOS","C9H9NO","C5H7NO","C3H5NO2","C4H7NO2",
               "C11H10N2O","C9H9NO2","C5H9NO","C4H4O2","C3H5O","C10H12N5O6P",
               "C10H15N2O3S","C10H14N2O2S","CH2ON","C21H34N7O16P3S",
               "C21H33N7O15P3S","C10H15N3O5S","C5H7","C3H2O3","C16H30O",
               "C8H8NO5P","CH3N2O","C5H4N5","C10H11N5O3","C10H13N5O9P2",
               "C10H12N5O6P","C9H13N3O10P2","C9H12N3O7P","C4H4N3O",
               "C10H13N5O10P2","C10H12N5O7P","C5H4N5O","C10H11N5O4",
               "C10H14N2O10P2","C10H12N2O4","C5H5N2O2","C10H13N2O7P",
               "C9H12N2O11P2","C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2",
               "C2H2O","C2H2","CO2","CHO2","H2O","H3O6P2","C2H4","CO","C2O2",
               "H2","O","P","C2H2O","CH2","HPO3","NH2","PP","NH","SO3","N",
               "C6H10O5","C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6",
               "C6H10O5","C18H30O15"], ncores=1):
    """
    Compute matrix of biochemical connections. Either based on a list of
    possible connections in the form of a list of formulas or based on the
    reactions present in the database.
    
    Parameters
    ----------
    DB: pandas dataframe containing the database against which the annotation
        is performed. The DB must contain the following columns in this exact
        order (optional fields can contain None):
        - id: unique id of the database entry (e.g., 'C00031') - necessary
        - name: compound name (e.g., 'D-Glucose') - necessary
        - formula: chemical formula (e.g., 'C6H12O6') - necessary
        - inchi: inchi string - optional
        - smiles: smiles string - optional
        - RT: if known, retention time range (in seconds) where this compound
                is expected to elute (e.g., '30;60') - optional
        - adductsPos: list of adducts that should be considered in positive mode
                      for this entry (e.g.,'M+Na;M+H;M+') - necessary
        - adductsNeg: list of adducts that should be considered in negative
                      mode for this entry (e.g.,'M-H;M-2H') - necessary
        - description: comments on the entry - optional
        - pk: previous knowledge on the likelihood of this compound to be
             present in the sample analyse. The value has to be between 1
             (compound likely to be present in the sample) and 0 (compound
             cannot be present in the sample).
        - MS2: id for the MS2 database entries related to this compound
               (optional)
        - reactions: list of reactions ids involving this compound
                    (e.g., 'R00010 R00015 R00028')-optional, but necessary if 
                    mode='reactions'.
    annotations: If equal to None (default) all entries in the DB are considered 
                (used to pre-compute the Bio matrix), alternatively it should be
                a dictionary containing all the possible annotations for the
                measured features. The keys of the dictionary are the unique ids
                for the features present in df. For each feature, the
                annotations are summarized in a pandas dataframe. Output of
                functions MS1annotation(), MS1annotation_Parallel(),
                MSMSannotation() or MSMSannotation_Parallel. In this case
                only the entries currently considered as possible annotations
                are used.
    mode: either 'reactions' (connections are computed based on the reactions
          present in the database) or 'connections' (connections are computed
          based on the list of connections provided). Default 'reactions'.
    connections: list of possible connections between compounds defined as
                formulas. Only necessary if mode='connections'. A list of
                common biotransformations is provided as default.
    ncores: default value 1. Number of cores used
    
    Returns
    -------
        Bio: dataframe containing all the possible connections computed.
    """
    if ncores==1:
        print("computing all possible biochemical connections")
        start = time.time()
        DB = DB.replace(numpy.nan,None)
        DB.loc[DB.reactions==None,'reactions']=''
            
        ### getting the ids of all possible hits to the database
        if annotations is None:
            all_ids = DB['id'].to_list()
        else:
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
                Bio = Bio+[iterations.bio_single_iter_connections(all_ids_DB,all_forms_DB,connections,x,y)]
    
    
        
        elif mode=='reactions':
            print("considering the reactions stored in the database ...")
            all_rs_DB = DB['reactions'].to_list()
            all_rs_DB= ['' if v is None else v for v in all_rs_DB]
            for x,y in itertools.combinations(all_ids, 2):
                Bio = Bio+[iterations.bio_single_iter_reactions(all_ids_DB,all_rs_DB,x,y)]
    
    
        
        else:
            print('ERROR: mode can only be connections or reactions')
            return
    
        Bio = [i for i in Bio if i != ('x','x')]
        Bio = pandas.DataFrame(Bio, index=None)
        end = time.time()   
        print(round(end - start,1), 'seconds elapsed')

    elif ncores>1:
        print("computing all possible biochemical connections - Parallelized")
        start = time.time()
        DB = DB.replace(numpy.nan,None)
        DB.loc[DB.reactions==None,'reactions']=''
        
        ### getting the ids of all possible hits to the database
        if annotations is None:
            all_ids = DB['id'].to_list()
        else:
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
            Bio = pool_obj.starmap(partial(iterations.bio_single_iter_connections,all_ids_DB,all_forms_DB,connections),itertools.combinations(all_ids, 2))
            pool_obj.terminate()
    
        
        elif mode=='reactions':
            print("considering the reactions stored in the database ...")
            all_rs_DB = DB['reactions'].to_list()
            all_rs_DB= ['' if v is None else v for v in all_rs_DB]
            pool_obj = multiprocessing.Pool(ncores)
            Bio = pool_obj.starmap(partial(iterations.bio_single_iter_reactions,all_ids_DB,all_rs_DB),itertools.combinations(all_ids, 2))
            pool_obj.terminate()
    
        
        else:
            print('ERROR: mode can only be connections or reactions')
            return
    
        
        Bio = [i for i in Bio if i != ('x','x')]
        Bio = pandas.DataFrame(Bio, index=None)
        end = time.time()   
        print(round(end - start,1), 'seconds elapsed')
    else:
        raise ValueError("ncores must be >=1")
        
    return(Bio)





def Gibbs_sampler_bio(df,annotations,Bio,noits=100,burn=None,delta_bio=1,
                      all_out=False,zs=None):
    """
    Gibbs sampler considering only possible biochemical connections. The
    function computes the posterior probabilities of the annotations
    considering the possible biochemical connections reported in Bio.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
    annotations: a dictionary containing all the possible annotations for the
                 measured features. The keys of the dictionary are the unique
                 ids for the features present in df. For each feature, the
                 annotations are summarized in a pandas dataframe. Output of
                 functions MS1annotation(), MS1annotation_Parallel(),
                 MSMSannotation() or MSMSannotation_Parallel
    Bio: dataframe (2 columns), reporting all the possible connections between
         compounds. It uses the unique ids from the database. It could be the
         output of Compute_Bio() or Compute_Bio_Parallel().
    noits: number of iterations if the Gibbs sampler to be run
    burn: number of iterations to be ignored when computing posterior
          probabilities. If None, is set to 10% of total iterations
    delta_bio: parameter used when computing the conditional priors.
               The parameter must be positive. The smaller the parameter the
               more weight the adducts connections have on the posterior
               probabilities. Default 1.
    all_out: logical value. If true the list of assignments found in each
            iteration is returned by the function. Default False.
    zs: list of assignments computed in a previous run of the Gibbs sampler.
        Optional, default None.
    
    Returns
    -------
    annotations: a dictionary containing all the possible annotations for the
                 measured features. The keys of the dictionary are the unique
                 ids for the features present in df. For each feature, the
                 annotations are summarized in a pandas dataframe.

    """
    df=df.replace('None',None)
    start = time.time()
    print("computing posterior probabilities including biochemical connections")
    print("initialising sampler ...")
    all_ids = []
    for k in annotations.keys():
        tmp = annotations[k]
        all_ids=all_ids+list(tmp['id'])
    all_ids=list(set(all_ids))

    ind = []
    for k in range(0,len(Bio.index)):
        if Bio.iloc[k,0] in all_ids and Bio.iloc[k,1] in all_ids:
            ind.append(k)
        
    Bio=Bio.iloc[ind,:]
    del all_ids
    del ind
    del tmp
    
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
        ca, ca_id = iterations.gibbs_sampler_bio_iter(indk,ks,annotations,Bio,ca_id,ca,delta_bio,it)
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






def Gibbs_sampler_bio_add(df,annotations,Bio,noits=100,burn=None,delta_bio=1,
                          delta_add=1,all_out=False,zs=None):
    """
    Gibbs sampler considering both biochemical and adducts connections. The
    function computes the posterior probabilities of the annotations
    considering the possible biochemical connections reported in Bio and the
    possible adducts connection.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
    annotations: a dictionary containing all the possible annotations for the
                 measured features. The keys of the dictionary are the unique
                 ids for the features present in df. For each feature, the
                 annotations are summarized in a pandas dataframe. Output of
                 functions MS1annotation(), MS1annotation_Parallel(),
                 MSMSannotation() or MSMSannotation_Parallel
    Bio: dataframe (2 columns), reporting all the possible connections between
         compounds. It uses the unique ids from the database. It could be the
         output of Compute_Bio() or Compute_Bio_Parallel().
    noits: number of iterations if the Gibbs sampler to be run
    burn: number of iterations to be ignored when computing posterior
          probabilities. If None, is set to 10% of total iterations
    delta_bio: parameter used when computing the conditional priors.
               The parameter must be positive. The smaller the parameter the
               more weight the adducts connections have on the posterior
               probabilities. Default 1.
    delta_add: parameter used when computing the conditional priors. The
               parameter must be positive. The smaller the parameter the more
               weight the adducts connections have on the posterior
               probabilities. Default 1.
    all_out: logical value. If true the list of assignments found in each
            iteration is returned by the function. Default False.
    zs: list of assignments computed in a previous run of the Gibbs sampler.
        Optional, default None.
    
    Returns
    -------
    annotations: the function modifies the annotations dictionary by adding 2
                columns to each entry. One named 'post Gibbs' contains the
                posterior probabilities computed. The other is called
                'chi-square pval' containing the p-value from a chi-squared
                test comparing the 'post' with the 'post Gibbs' probabilities.
    zs: optional, if all_out==True, the function return the full list of
        assignments computed. This allows restarting the sampler from where you
        are from a previous run
    """
    start = time.time()
    print("computing posterior probabilities including biochemical and adducts connections")
    print("initialising sampler ...")
    df=df.replace('None',None)
    all_ids = []
    for k in annotations.keys():
        tmp = annotations[k]
        all_ids=all_ids+list(tmp['id'])
    all_ids=list(set(all_ids))

    ind = []
    for k in range(0,len(Bio.index)):
        if Bio.iloc[k,0] in all_ids and Bio.iloc[k,1] in all_ids:
            ind.append(k)
        
    Bio=Bio.iloc[ind,:]
    del all_ids
    del ind
    del tmp
 
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
        ca, ca_id = iterations.gibbs_sampler_bio_add_iter(indk,ks,rids,annotations,Bio,ca_id,ca,delta_bio,delta_add,it)
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


def simpleIPA(df,ionisation,DB,adductsAll,ppm,dfMS2=None,DBMS2=None,noits=100,
              burn=None,delta_add=None,delta_bio=None,Bio=None,
              mode='reactions',CSunk=0.5,isodiff=1,ppmiso=100,ncores=1,
              me=5.48579909065e-04,ratiosd=0.9,ppmunk=None,ratiounk=None,
              ppmthr=None,pRTNone=None,pRTout=None,mzdCS=0, ppmCS=10,
              evfilt=False,
              connections = ["C3H5NO", "C6H12N4O", "C4H6N2O2", "C4H5NO3",
                             "C3H5NOS", "C6H10N2O3S2","C5H7NO3","C5H8N2O2",
                             "C2H3NO","C6H7N3O","C6H11NO","C6H11NO","C6H12N2O",
                             "C5H9NOS","C9H9NO","C5H7NO","C3H5NO2","C4H7NO2",
                             "C11H10N2O","C9H9NO2","C5H9NO","C4H4O2","C3H5O",
                             "C10H12N5O6P","C10H15N2O3S","C10H14N2O2S","CH2ON",
                             "C21H34N7O16P3S","C21H33N7O15P3S","C10H15N3O5S",
                             "C5H7","C3H2O3","C16H30O","C8H8NO5P","CH3N2O",
                             "C5H4N5","C10H11N5O3","C10H13N5O9P2",
                             "C10H12N5O6P","C9H13N3O10P2","C9H12N3O7P",
                             "C4H4N3O","C10H13N5O10P2","C10H12N5O7P","C5H4N5O",
                             "C10H11N5O4","C10H14N2O10P2","C10H12N2O4",
                             "C5H5N2O2","C10H13N2O7P","C9H12N2O11P2",
                             "C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2",
                             "C2H2O","C2H2","CO2","CHO2","H2O","H3O6P2","C2H4",
                             "CO","C2O2","H2","O","P","C2H2O","CH2","HPO3",
                             "NH2","PP","NH","SO3","N","C6H10O5",
               "C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6","C6H10O5",
               "C18H30O15"]):
    """
    Wrapper function performing the whole IPA pipeline.
    
    Parameters
    ----------
    df: pandas dataframe containing the MS1 data. It should be the output of the
        function ipa.map_isotope_patterns()
        
    DB: pandas dataframe containing the database against which the annotation
        is performed. The DB must contain the following columns in this exact
        order (optional fields can contain None):
        - id: unique id of the database entry (e.g., 'C00031') - necessary
        - name: compound name (e.g., 'D-Glucose') - necessary
        - formula: chemical formula (e.g., 'C6H12O6') - necessary
        - inchi: inchi string - optional
        - smiles: smiles string - optional
        - RT: if known, retention time range (in seconds) where this compound
                is expected to elute (e.g., '30;60') - optional
        - adductsPos: list of adducts that should be considered in positive mode
                      for this entry (e.g.,'M+Na;M+H;M+') - necessary
        - adductsNeg: list of adducts that should be considered in negative
                      mode for this entry (e.g.,'M-H;M-2H') - necessary
        - description: comments on the entry - optional
        - pk: previous knowledge on the likelihood of this compound to be
             present in the sample analyse. The value has to be between 1
             (compound likely to be present in the sample) and 0 (compound
             cannot be present in the sample).
        - MS2: id for the MS2 database entries related to this compound
               (optional)
        - reactions: list of reactions ids involving this compound
                    (e.g., 'R00010 R00015 R00028')-optional, but necessary if 
                    mode='reactions'.
    adductsAll:a dataframe containing information on all possible adducts. 	        1
    ppm: accuracy of the MS instrument used
    dfMS2: pandas dataframe containing the MS2 data (optional). It must contain
           3 columns:
                   -id: an unique id for each feature for which the MS2 spectrum
                       was acquired (same as in df)
                   -spectrum: string containing the spectrum inforamtion in the
                              following format 'mz1:Int1 mz2:Int2 mz3:Int3 ...'
                   -ev: collision energy used to aquire the fragmentation
                       spectrum
    DBMS2: pandas dataframe containing the database containing the MS2
           information (optional)
    evfilt: Default value False. If true, only spectra acquired with the same
            collision energy are considered.
    noits: number of iterations if the Gibbs sampler to be run
    burn: number of iterations to be ignored when computing posterior
          probabilities. If None, is set to 10% of total iterations
    delta_bio: parameter used when computing the conditional priors.
               The parameter must be positive. The smaller the parameter the
               more weight the adducts connections have on the posterior
               probabilities. Default 1.
    delta_add: parameter used when computing the conditional priors. The
               parameter must be positive. The smaller the parameter the more
               weight the adducts connections have on the posterior
               probabilities. Default 1.
    Bio: dataframe (2 columns), reporting all the possible connections between
         compounds. It uses the unique ids from the database. It could be the
         output of Compute_Bio() or Compute_Bio_Parallel().
    mode: either 'reactions' (connections are computed based on the reactions
          present in the database) or 'connections' (connections are computed
          based on the list of connections provided). Default 'reactions'.
   CSunk: cosine similarity score associated with the 'unknown' annotation.
           Default 0.7    
   isoDiff: Default value 1. Difference between isotopes of charge 1, does not
            need to be exact
    ppmiso: Default value 100. Maximum ppm value allowed between 2 isotopes.
            It is very high on purpose
    ncores: default value 1. Number of cores used
    

    me: accurate mass of the electron. Default 5.48579909065e-04
    ratiosd: default 0.9. It represents the acceptable ratio between predicted
            intensity and observed intensity of isotopes. it is used to compute
            the shape parameters of the lognormal distribution used to
            calculate the isotope pattern scores as sqrt(1/ratiosd)
    ppmunk: ppm associated to the 'unknown' annotation. If not provided equal
            to ppm.
    ratiounk: isotope ratio associated to the 'unknown' annotation. If not
              provided equal to 0.5
    ppmthr: Maximum ppm possible for the annotations. Ff not provided equal to
            2*ppm
    pRTNone: Multiplicative factor for the RT if no RTrange present in the
            database. If not provided equal to 0.8
    pRTout: Multiplicative factor for the RT if measured RT is outside the
            RTrange present in the database. If not provided equal to 0.4
    mzdCS: maximum mz difference allowed when computing cosine similarity
           scores. If one wants to use this parameter instead of ppmCS, this
           must be set to 0. Default 0.
    ppmCS: maximum ppm allowed when computing cosine similarity scores.
           If one wants to use this parameter instead of mzdCS, this must be
           set to 0. Default 10.
    connections: list of possible connections between compounds defined as
                formulas. Only necessary if mode='connections'. A list of
                common biotransformations is provided as default.
    Output:
        annotations: a dictionary containing all the possible annotations for the measured features. The keys of the dictionary are the
                     unique ids for the features present in df. For each feature, the annotations are summarized in a pandas dataframe.
    """
    df=df.replace('None',None)
    # mapping isotopes
    if len(df.columns)==5:
        map_isotope_patterns(df,isoDiff=isodiff, ppm=ppmiso,ionisation=ionisation)
    elif len(df.columns)==8:
        print("isotopes already mapped")
    else:
        raise ValueError("df not in the correct format")
        
    
    # computing all adducts
   
    allAdds = compute_all_adducts(adductsAll= adductsAll, DB=DB, ionisation=ionisation ,ncores=ncores)
   

    # computing priors
    if (dfMS2 is None) or (DBMS2 is None):
        annotations = MS1annotation(df=df,allAdds=allAdds,ppm=ppm,me=me,ratiosd=ratiosd,ppmunk=ppmunk,
        ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout,ncores=ncores)
    else:
    
        annotations = MSMSannotation(df=df,dfMS2=dfMS2,allAdds=allAdds,DBMS2=DBMS2,ppm=ppm,me=me,ratiosd=ratiosd,
        ppmunk=ppmunk,ratiounk=ratiounk,ppmthr=ppmthr,pRTNone=pRTNone,pRTout=pRTout,mzdCS=mzdCS,ppmCS=ppmCS,
        CSunk=CSunk,evfilt=evfilt,ncores=ncores)

    # computing Bio matrix (if necessary)
    if (Bio is None) and (delta_bio is not None):
        Bio=Compute_Bio(DB=DB,annotations=annotations,mode=mode,connections=connections,ncores=ncores)
        
    # Gibbs sampler (if needed). Which one based on the inputs
    if (Bio is not None) and (delta_bio is not None) and (delta_add is not None):
        Gibbs_sampler_bio_add(df=df,annotations=annotations,Bio=Bio,noits=noits,burn=burn,delta_bio=delta_bio,delta_add=delta_add)
    elif (Bio is not None) and (delta_bio is not None) and (delta_add is None):
        Gibbs_sampler_bio(df=df,annotations=annotations,Bio=Bio,noits=noits,burn=burn,delta_bio=delta_bio)
    elif (Bio is None) and (delta_bio is None) and (delta_add is not None):
        Gibbs_sampler_add(df=df,annotations=annotations,noits=noits,burn=burn,delta_add=delta_add)
    
    return(annotations)

        

