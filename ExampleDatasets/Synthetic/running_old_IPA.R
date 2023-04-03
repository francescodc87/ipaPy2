setwd('~/Desktop/ipaPy2-main/ExampleDatasets/Synthetic/')
library(IPA)
data("isotopes")
dfPOS <- read.csv('positive_synth_dataset.csv', stringsAsFactors = F)
load('DB_IPA1.RData')

ionisation="positive"
ppm=5

Hits <- find.hits(adducts.matrix= all_adducts_POS,
                  dataset=dfPOS[,3:5], ppm.thr= 30,
                  RTwin=60,relation.id = dfPOS$rel.id,
                  isotopes=isotopes, 
                  iso.threshold=1)


Prior <- compute.Priors(Hits=Hits, dataset=dfPOS[,3:5],
                        pk=rep(1,nrow(Hits$all.formulas)),
                        ppm=ppm,unknown.ppm = 5,pr.lim = 0)


### building ADD matrix
ADD <- build.add.connenctivity.matrix(Prior=Prior,  DB=DB,
                                      ionisation,fully.connected=TRUE)

### building ISO matrix
ISO <- build.iso.connectivity.matrix(Prior=Prior, DB=DB, ratios=TRUE)

### building BIO matrix
BIO1 <- build.bio.connectivity.matrix(Prior=Prior, DB=DB,
                                      ionisation,
                                      connection.type="reactions")


Post <- IPAposteriors(P=Prior,Iso = ISO, Add = ADD, Bio = BIO1,
                               Int = as.numeric(dfPOS[,3]), ratio.toll = 0.8,
                               delta.iso = .1, delta.add = .1, delta.bio = .1,
                               allsamp = T, no.its = 5000, burn = 1000,
                               rel.id =dfPOS$rel.ids)


Final.res <- ParseIPAresults(Post=Post, Prior=Prior,dataset =dfPOS[,3:5], DB = DB,
                             IDs=dfPOS$ids)

ResultsPOS <- list(Hits=Hits, Prior= Prior, Post=Post, Final.res=Final.res)


postIPA1 <- rep(0,nrow(dfPOS))
for(k in 1:nrow(dfPOS)){
  id <- dfPOS$ids[k]
  comp <- dfPOS$Compound.ID[k]
  form <- dfPOS$Formula[k]
  tmp <- ResultsPOS$Final.res[[id]]
  ind <- which(tmp[,1]==comp & tmp[,4]==form)
  if(length(ind)>1){
    cat(k)
  }
  postIPA1[k] <- as.numeric(tmp[ind,9])
}

dfPOS$'post IPA1' <- postIPA1
ind <- which(dfPOS$isotope=='mono')
sum(log(dfPOS$`post IPA1`[ind]))
write.csv(dfPOS[ind,c(1:8,11)], 'res_IPAV1_POS.csv', row.names = F)





dfNEG <- read.csv('negative_synth_dataset.csv', stringsAsFactors = F)

ionisation="negative"
ppm=5

Hits <- find.hits(adducts.matrix= all_adducts_NEG,
                  dataset=dfNEG[,3:5], ppm.thr= 30,
                  RTwin=60,relation.id = dfNEG$rel.id,
                  isotopes=isotopes, 
                  iso.threshold=1)


Prior <- compute.Priors(Hits=Hits, dataset=dfNEG[,3:5],
                        pk=rep(1,nrow(Hits$all.formulas)),
                        ppm=ppm,unknown.ppm = 5,pr.lim = 0)


### building ADD matrix
ADD <- build.add.connenctivity.matrix(Prior=Prior,  DB=DB,
                                      ionisation,fully.connected=TRUE)

### building ISO matrix
ISO <- build.iso.connectivity.matrix(Prior=Prior, DB=DB, ratios=TRUE)

### building BIO matrix
BIO1 <- build.bio.connectivity.matrix(Prior=Prior, DB=DB,
                                       ionisation,
                                       connection.type="reactions")


Post <- IPAposteriors(P=Prior,Iso = ISO, Add = ADD, Bio = BIO1,
                      Int = as.numeric(dfNEG[,3]), ratio.toll = 0.8,
                      delta.iso = .1, delta.add = .1, delta.bio = .1,
                      allsamp = T, no.its = 5000, burn = 1000,
                      rel.id =dfNEG$rel.ids)


Final.res <- ParseIPAresults(Post=Post, Prior=Prior,dataset =dfNEG[,3:5], DB = DB,
                             IDs=dfNEG$ids)
ResultsNEG <- list(Hits=Hits, Prior= Prior, Post=Post, Final.res=Final.res)





postIPA1 <- rep(0,nrow(dfNEG))
for(k in 1:nrow(dfNEG)){
  id <- dfNEG$ids[k]
  comp <- dfNEG$Compound.ID[k]
  form <- dfNEG$Formula[k]
  tmp <- ResultsNEG$Final.res[[id]]
  ind <- which(tmp[,1]==comp & tmp[,4]==form)
  if(length(ind)==0){
    cat(k)
  }
  postIPA1[k] <- as.numeric(tmp[ind,9])
}

dfNEG$'post IPA1' <- postIPA1
ind <- which(dfNEG$isotope=='mono')
sum(log(dfNEG$`post IPA1`[ind]))
write.csv(dfNEG[ind,c(1:8,11)], 'res_IPAV1_NEG.csv',row.names = F)




