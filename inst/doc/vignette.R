source("funcsToUpdate_1_2017.R")
source("old_function.R")

load("AllExamples.RData")
#Vaccine data used in CellCODE paper
allPaths=combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,  canonicalPathways)
#for larger datasets litvakSigs may also be useful
#allPaths=combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,  canonicalPathways, litvakSigs)
cm.genes=commonRows(allPaths, vacData)
vacDataN=rowNorm(vacData)
num.pc(vacDataN[cm.genes,])
#the result is 24
plierResult=PLIER(vacDataN[cm.genes,], allPaths[cm.genes,],k=24, trace=T, max.iter=150)
#Correlate with SPVs from cellCODE
plotMat(cor(t(plierResult$B), SPVs))
#we have nice one-to-one  correspondence, though the "DendriticCell" signature from CellCODE is more closely related to the Type-I interferon transcriptional response so it is probably not cell-type induced variation.
#Visualize the cross-validation results
plotMat(plierResult$Uauc)

#plot all of U 
plotU(plierResult,auc.cutoff = 0.5, pval.cutoff = 1)
#visualize the top genes
plotTopZ(plierResult, vacDataN, allPaths, top = 10)

#the "PID_ATF2_PATHWAY" looks a little tenuous
##we can check it's stats
plierResult$summary[which(plierResult$summary$`LV index`==3),]
#the association with "PID_ATF2_PATHWAY" is not significant;  this pathway has only 42 genes that are also present in the dataset which makes achieving significance harder 
#In general, the signaling pathway  doesn't have perfect correspondence to the transcriptional signature, but the transcriptional signature can still be of interest and we can plot it by itself with more genes
plotTopZ(plierResult, vacDataN, allPaths, index=c(3), top=50)
#looks like an AP1 pathway that regulates IL8 transcription possibly downstream of TNF signaling (we have FOS, JUNB, FOSB (not actually a memeber of the PID_ATF2_PATHWAY geneset))


##HCC example
##
CancerPath=combinePaths(canonicalPathways, chemgenPathways, oncogenicPathways)
cmHCC=commonRows(HCCdataTumor, CancerPath)
#remove small pathways, not strictly necessary but saves computation time 
#by making the pathway/geneset matrix smaller
ii=which(colSums(CancerPath[cmHCC,])<20)
HCCpath=CancerPath[, -ii]
#prescale the data
HCCdataN=rowNorm(HCCdataTumor[cmHCC,])

#Precompute Chat, this is used to define active pathways and is expensive for large pathway sets
#helpful if we want to run plier with different parameters
HCCchat=computeChat(CancerPath[cmHCC,])
#num.pc is also expensive, for this dataset it is 52
plierResultHCC=PLIER(HCCdataN, CancerPath[cmHCC,], k = 52, Chat = HCCchat, trace=T)
#plot with a high AUC cutoff so it is not too busy
plotU(plierResultHCC, auc.cutoff = 0.9)
#we found two immune components, interferon alpha and genes related to interferon gamma/CD8/Th1 response
plotTopZ(plierResultHCC, HCCdataN, CancerPath, index=c(26, 40), top = 20)
#note component 40 has many pathways, it is named with the top pathway but "inPathway" refers to their union




