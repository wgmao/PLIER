library(PLIER)
data(HCCdata)
data(canonicalPathways)
data(chemgenPathways)
data(oncogenicPathways)


HCCpath=combineGSmats(canonicalPathways, chemgenPathways, oncogenicPathways)
cmHCC=commonRows(HCCdata, HCCpath)

#remove small pathways, not strictly necessary but saves computation time 
#by making the pathway/geneset matrix smaller
ii=which(colSums(HCCpath[cmHCC,])<20)
HCCpath=HCCpath[, -ii]

dim(HCCpath)
#HCCdata how genes in rows and samples in columns
dim(HCCdata)
#use only the samples with survival info: i.e. the tumor samples

#also presecale the data
HCCdataUse=rowNorm(HCCdata[cmHCC,cmSurvival])
#precompute the SVD, useful if testing different parameters
svdresHCC=svd(HCCdataUse)


#compute the number of components
#num.pc(HCCdataUse)
#the result if 52
#setting max.iter=250, L1=25, k=52, all other parameters are default
plierRes=PLIER(data=HCCdataUse, HCCpath[cmHCC,],svdres=svdresHCC, k=52, L1=25, max.iter=250, trace=T)
plotMat(plierRes$U)
#check for a CTNNB1 mutation correlation
#these are activating mutations so excpect a positive correlation
data(CTNNB1mut)
plot(sort(corRes<-cor(t(plierRes$B), CTNNB1mut[cmSurvival], use = "p")), xlab="LV", ylab="Correlation with CTNNB1 mutation")
#and the correlated LV is
iBest=which(corRes==max(corRes))
#check the associated pathways
which(plierRes$U[,iBest]>0)

#check survival
survRaw=testSurvival(plierRes$B, survTime, survStatus)
#include the first 10 PCs from the SVD to correct for global subtype structure
survCor=testSurvival(plierRes$B, survTime, survStatus, t(svdresHCC$v[, 1:10]))

#plot survival predictors
iiSurv=which(survRaw$q.value<0.2)
plotMat(plierRes$U[,iiSurv])


#plot structure corrected survival predictors
iiSurv=which(survCor$q.value<0.2)
plotMat(plierRes$U[,iiSurv])
