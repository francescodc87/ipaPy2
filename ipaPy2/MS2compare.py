### functions necessary for the comparison of fragmentation spectra
import pandas
import numpy
import math
import statistics

__author__ = "Francesco Del Carratore, Juraj Borka"
__maintainer__ = "Francesco Del Carratore"
__email__ = "francesco.delcarratore@gmail.com"

class SampleFragmentsItem():

    def __init__(self,mz,intensity):
        self.mz = mz
        self.intensity = intensity

    def length(self):
        return(len(self.mz))


def diff(x):
    counter = 0
    out = []
    while counter < len(x)-1:
        out.append(x[counter+1]-x[counter])
        counter += 1
    return out

def cumsum(x):
    out = []
    sum = 0
    for i in x:
        sum += i
        out.append(sum)
    return out


def group_mz_values(x,mzd,ppm = 10):
    mzdiff = diff(x)
    input_array = [0]

    counter = 0
    if int(ppm) > 0:
        while counter < len(mzdiff):
            if mzdiff[counter] >= (mzd + x[-len(x)] * int(ppm) / 1e6):
                input_array.append(1)
            else:
                input_array.append(0)
            counter += 1
        input_array.append(1)
        return cumsum(input_array)
    else:
        while counter < len(mzdiff):
            if mzdiff[counter] >= mzd:
                input_array.append(1)
            else:
                input_array.append(0)
            counter += 1
        input_array.append(1)
        return cumsum(input_array)


def split(x,mapping):
    outDict = {}
    counter = 0
    while counter < len(x):
        if mapping[counter] in outDict:
            outDict[mapping[counter]].append(x[counter])
        else:
            outDict[mapping[counter]] = ([x[counter]])
        counter += 1

    return outDict


def lengths(x):
    out = []
    for key in x:
        out.append((key,len(x[key])))
    return out

def ConsensusSpec(x,minProp,mzd = 0,ppm = 10):
    if len(x) == 1:
        highest = max(x[0].intensity)
        new_ints = []
        for j in x[0].intensity:
            new_ints.append((j/highest)*100)
        x[0].intensity = new_ints
        return(x[0])

    allIntensities = []

    #normalise Data

    for i in x:
        highest = max(i.intensity)
        new_ints = []
        for j in i.intensity:
            new_ints.append((j/highest)*100)
        i.intensity = new_ints

    for i in x:
        allIntensities += i.intensity

    mzList = []
    for sample in x:
        mzList.append(sample.mz)

    mzListLengths = []
    for mzs in mzList:
        mzListLengths.append(len(mzs))

    tracker = []
    counter = 1
    for i in mzList:
        for j in i:
            tracker.append((counter,j,allIntensities[counter-1]))
            counter += 1
    tracker.sort(key=lambda tup: tup[1])

    allMzOrder = []
    allMz = []
    intensities = []
    for i in tracker:
        allMzOrder.append(i[0])
        allMz.append(i[1])
        intensities.append(i[2])

    mz_groups = group_mz_values(allMz,mzd,ppm)

    mz_split = split(allMz,mz_groups)
    int_split = split(intensities,mz_groups)

    keep = []
    for i in mz_split:
        if len(mz_split[i]) >= len(x)*minProp:
            keep.append((i,True))
        else:
            keep.append((i,False))


    outMZS = []
    outINTS = []

    allCounter = 0
    while allCounter < len(mz_split):
        if keep[allCounter][1] == True:
            outMZS.append(sum(mz_split[allCounter])/len(mz_split[allCounter]))
            outINTS.append(sum(int_split[allCounter])/len(int_split[allCounter]))
        allCounter += 1

    #Finish normalise

    out = SampleFragmentsItem(outMZS,outINTS)

    return out

def compareSpectra(a,b):
    if (len(a.mz) == len(a.intensity)) and (len(b.mz) == len(b.intensity)):
        
        binnedSpectra = binSpectra(a,b)
        inten = {}
        counter = 1
        for i in binnedSpectra:
            inten[str(counter)] = i["intensity"]
            counter += 1

        out = dotProduct(inten["1"],inten["2"])
        return out
    else:
        return "unmatched mz/intensity"


def dotProduct(x,y):
    mat = numpy.matmul(x,y)

    sqrtx = 0
    for i in x:
        sqrtx += i*i

    sqrty = 0
    for i in y:
        sqrty += i*i

    sqrtx = math.sqrt(sqrtx)
    sqrty = math.sqrt(sqrty)

    return mat/(sqrtx*sqrty)

def binSpectra(x,y):
    holder = x.mz + y.mz
    floor = math.floor(min(holder))
    ceiling = math.ceil(max(holder))

    breaks = list(range(floor,ceiling+1))
    breaks = fixBreaks(breaks,[min(holder),max(holder)])

    return [binSpectrum(x,breaks),binSpectrum(y,breaks)]


def binSpectrum(x,breaks):
    bins = binValues(x.intensity, x.mz, breaks)

    obj = {}
    obj["mz"] = bins["mid"]
    obj["intensity"] = bins["x"]
    obj["tic"] = sum(obj["intensity"])
    obj["peaksCount"] = len(obj["mz"])

    return obj


def binValues(ints,mzs,brks):
    breaks = fixBreaks(brks,[min(mzs),max(mzs)])
    nbrks = len(breaks)
    idx = findInterval(mzs,breaks)
    iints = []
    for i in range(1,nbrks):
        iints.append(0)
    s = split(ints,idx)
    for j in s:
        s[j] = max(s[j])
    for i in s:
        iints[i-1] = s[i]

    out = {"x":[],"mid":[]}
    x = []
    mid = []

    counter = 0

    while counter < len(iints):
        try:
            x.append(iints[counter])
            mid.append(statistics.mean([breaks[counter],breaks[counter+1]]))
        except:
            pass
        counter += 1

    out["x"] = x
    out["mid"] = mid

    return out




def fixBreaks(brks, rng):
    if brks[len(brks)-1] <= rng[1]:
        out = brks
        out.append(max([rng[1]+1e-6,brks[len(brks)-1] + statistics.mean(diff(brks))]))
        return out
    else:
        return brks


def findInterval(x,vec):
    vec.append(max(x)+1)
    tracker = 0
    out = []
    for i in x:
        try:
            if i >= vec[tracker]:
                tracker += 1
                out.append(tracker)
            else:
                out.append(tracker)
        except:
            pass
    vec.pop()
    return out


def cosine_similarity(DBsp,Msp,mzd=0,ppm=10):
    ## extract data from string
    DBsp = DBsp.split()
    Msp = Msp.split()
    mz1 = [0]*len(DBsp)
    int1 = [0]*len(DBsp)
    for k in range(0,len(mz1)):
        tmp = DBsp[k].split(':')
        mz1[k] = float(tmp[0])
        int1[k]=float(tmp[1])
    mz2 = [0]*len(Msp)
    int2 = [0]*len(Msp)
    for k in range(0,len(mz2)):
        tmp = Msp[k].split(':')
        mz2[k]=float(tmp[0])
        int2[k]=float(tmp[1])
    
    
    allMz = mz1+mz2
    allInts = int1+int2
    splab = ['A']*len(mz1)+['B']*len(mz2)
    spdf = pandas.DataFrame([allMz,allInts,splab]).transpose()
    spdf = spdf.sort_values(by=[0], ascending=True)

    mz_groups = group_mz_values(list(spdf[0]),mzd=mzd, ppm=ppm)
    mz_split = split(list(spdf[0]),mz_groups)
    lab_split = split(list(spdf[2]),mz_groups)
    int_split = split(list(spdf[1]),mz_groups)

    ks = mz_split.keys()

    mz_new = [0]*len(mz_split)
    Int_new1 = [0]*len(mz_split)
    Int_new2 = [0]*len(mz_split)

    counter = 0
    for k in ks:
        mz_new[counter] = sum(mz_split[k])/len(mz_split[k])
        intA = [int_split[k][i] for i,val in enumerate(lab_split[k]) if val=='A']
        intB = [int_split[k][i] for i,val in enumerate(lab_split[k]) if val=='B']
        Int_new1[counter] = sum(intA)
        Int_new2[counter] = sum(intB)
        counter=counter+1  
    Int_new1 = (numpy.array(Int_new1)/max(Int_new1))*100
    Int_new2 = (numpy.array(Int_new2)/max(Int_new2))*100
    sp1 = SampleFragmentsItem(numpy.array(mz_new),Int_new1)
    sp2 = SampleFragmentsItem(numpy.array(mz_new),Int_new2)
    out = compareSpectra(sp1,sp2)
    return(out)

