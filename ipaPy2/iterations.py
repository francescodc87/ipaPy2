from __future__ import annotations
import pandas
import molmass
from scipy import stats
import math
import random
from ipaPy2 import MS2compare
from ipaPy2 import util

__author__ = "Francesco Del Carratore"
__maintainer__ = "Francesco Del Carratore"
__email__ = "francesco.delcarratore@gmail.com"

def all_adducts_iter(DB,adductsAll,ionisation,db):
    f = DB['formula'][db]
    f = molmass.Formula(f)
    M=f.isotope.mass
    recs = []
    if ionisation==1:
        adds = DB['adductsPos'][db]
    else:
        adds = DB['adductsNeg'][db]
    adds = adds.split(';')
    adducts = adductsAll[adductsAll['name'].isin(adds)]
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


def MSMS_ann_iter1(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k):
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
                    DBsp = DBMS2[(DBMS2['compound_id']==MS2id) & (DBMS2['precursorType']==precType) & (DBMS2['collision.energy']==ev)]['spectrum']
                    #DBsp = list(DBMS2[(DBMS2['compound_id']==MS2id) & (DBMS2['precursorType']==precType)]['spectrum'])
                    if len(DBsp)>0:
                        CStmp =[]
                        for spDB in DBsp:
                            CStmp.append(MS2compare.cosine_similarity(spDB,Msp,mzdCS,ppmCS))
                        
                        CS.append(max(CStmp))
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



def MSMS_ann_iter2(df,dfMS2,allAdds,DBMS2,ppm,me,ratiosd,ppmthr,ppmunk,ratiounk,pRTNone,pRTout,mzdCS,ppmCS,CSunk,sigmaln,k):
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
                    #ev = Msps.iloc[s,2]
                    precType = allAdds['adduct'][hits[h]]
                    #DBsp = DBMS2[(DBMS2['compound_id']==MS2id) & (DBMS2['precursorType']==precType) & (DBMS2['collision.energy']==ev)]['spectrum']
                    DBsp = list(DBMS2[(DBMS2['compound_id']==MS2id) & (DBMS2['precursorType']==precType)]['spectrum'])
                    if len(DBsp)>0:
                        CStmp =[]
                        for spDB in DBsp:
                            CStmp.append(MS2compare.cosine_similarity(spDB,Msp,mzdCS,ppmCS))
                        
                        CS.append(max(CStmp))
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
            if idcp!='Unknown':
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
            if idcp!='Unknown':
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
