require(RColorBrewer)
require(gplots)
require(pheatmap)
require(glmnet)
require(rsvd)
require(qvalue)


#' @keywords  internal
#' Solves for the U coefficients making efficient utilizatoin of the lasso path
#' @param Z current Z estimate
#' @param Chat the inverse of the C matrix
#' @param priorMat the prior pathway or C matrix
#' @param penalty.factor Penalties for different pathways, must have size ncol(priorMat).  
#' @param pathwaySelection Method to use for pathway selection. 
#' @param glm_alpha The elsatic net alpha parameter
#' @param maxPath The maximum number of pathways to consider
#' @param target.frac The target fraction on non-zero columns of
#' @param L3 Solve with a given L3, no search
solveU=function(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, target.frac=0.7, L3=NULL){
  
  
  Ur=Chat%*%Z #get U by OLS
  Ur=apply(-Ur,2,rank) #rank
  Urm=apply(Ur,1,min)
  
  U=matrix(0,nrow=ncol(priorMat), ncol=ncol(Z))
  if(is.null(L3)){

    lambdas=exp(seq(-4,-12,-0.125))
    results=list()
    lMat=matrix(nrow=length(lambdas), ncol=ncol(Z))
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet(y=Z[,i], x=priorMat[,iip], penalty.factor = penalty.factor[iip], alpha=glm_alpha, lower.limits=0, lambda = lambdas,intercept=T,  standardize=F )
      #plot(gres)
      gres$iip=iip
      lMat[,i]=colSums(as.matrix(gres$beta)>0)
      results[[i]]=gres
    }
    fracs=rowMeans(lMat>0)
    iibest=which.min(abs(target.frac-fracs))
    iibest
    
    
    for(i in 1:ncol(Z)){
      U[results[[i]]$iip,i]=results[[i]]$beta[,iibest]
    }#for i
    rownames(U)=colnames(priorMat)
    colnames(U)=1:ncol(Z)
  
 Utmp=solveU(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10,  L3=lambdas[iibest])
 
 #stop()
    return(list(U=U, L3=lambdas[iibest]))  
  }
  else{ #do one fit with a given lambda
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet(y=Z[,i], x=priorMat[,iip], penalty.factor = penalty.factor[iip], alpha=glm_alpha, lower.limits=0, lambda = L3,intercept=T,  standardize=F )
      U[iip,i]=as.numeric(gres$beta)
    }
   
    return(U)
  }
}


#' @keywords  internal
wrapString=function(string, width=30){
  
string=lapply(string, function(s){ paste(strwrap(gsub("_", " ",s), width=width), collapse="\n")})
unlist(string)
}
#' @keywords  internal
QV=function(pval){
  
  x=try(qvalue(pval))
  
  if(!is.list(x)){
    warning("Q-value error, defaulting to BH")
    #hist(pval)
    return(p.adjust(pval, method="BH"))
  }
  else{
    return(x$qvalue)
  }
}

#' @keywords  internal
mydist = function(x) {
  as.dist(1 - t(cor(t(x))))
}

#' @keywords  internal
BH= function(pval){p.adjust(pval, method="BH")}


#' SVD based smoothing for single cell RNAseq data
#' @param svdres svd result 
#' @param k number of components to use
#' @export
DataSmooth=function(svdres,k){
    k=1:k
  ds=sweep(svdres$u[, k],2,svdres$d[k],"*")%*%t(svdres$v[, k])
  ds
}
#' Rename pathway matrix gene names. Useful for species conversion
#' @param pathway pathway matrix
#' @param map Gene name map. A single column data.frame or matrix with map-from in row names and map-to in the first column
#' @export
mapPathway=function(pathway, map){
  cm=commonRows(map, pathway)
  show(length(cm))
  pathway=pathway[cm,]
  rownames(pathway)=map[cm,1]
  pathway
}

#' Creates a binary cell-type marker matrix using prior results. This matrix can be used for other downstream tasks that require cell-type markers, such as CellCODE
#' @param plierRes A PLIER result
#' @param priorMat the binary prior information matrix that was used to compute the plierResult. Including this insures that only the genes annotated to the pathway(s) are included
#' @param num The number of marker genes to produce
#' @param index The indecies of PLIER latent variables that are believed to represent cell-type proportions (as opposed to other sources of correlation)
#' @export
plierResToMarkers=function(plierRes, priorMat, num=20, index=NULL){
  
  ii=which(colSums(plierRes$U)>0)
  if(! is.null(index)){
    ii=intersect(ii,index)
  }
  Zuse=plierRes$Z[,ii, drop=F]

  for(i in 1:length(ii)){
    lv=ii[i]
paths=names(which(plierRes$U[,lv]<0.01))
genes=names(which(rowSums(priorMat[,paths])>0))
genesNotInPath=setdiff(rownames(Zuse), genes)
Zuse[genesNotInPath,i]=0
}


    
    tag=apply(-Zuse,2,rank)
    colnames(tag)=rownames(plierRes$B)[ii]  
  iim=apply(tag,1,min)
  iig=which(iim<=num)
  tag=tag[iig,]
  iin=rowSums(tag<=num)
  iimulti=which(iin>1)
 if(length(iimulti)>0){
  message(paste0("Genes not matched uniquely: ", paste(names(iimulti), collapse=", ")))
}
  tag=(tag<=num)+1-1
  
  tag
  
}

#' combines binary pathway matricies into one, rownames are matched by name
#' 
#' @param ... binary pathway matricies
#' @export
combinePaths=function(...){
  cats=list(...)
  genes=character()
  
  for( i in 1:length(cats)){
    genes=c(genes, rownames(cats[[i]]))
    cats[[i]]=as.data.frame(cats[[i]])
  }
  
  genes=unique(genes)
  #show(genes)
  mat=matrix(nrow=length(genes), ncol=0)
  rownames(mat)=genes
  for( i in 1:length(cats)){
    cats[[i]]=as.matrix(cats[[i]][genes,])
    mat=cbind(mat, cats[[i]])
  }
  mat[is.na(mat)]=0
  
  mat
}

#' @keywords internal
AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  myres=list()
  if(posn>0&negn>0){
    res=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE);
   
    myres$low=res$conf.int[1]
    myres$high=res$conf.int[2]
    myres$auc=(res$statistic)/(posn*negn)
    myres$pval=res$p.value
  }
  else{
    myres$auc=0.5
    myres$pval=NA
  }
  return(myres)
}
#' get the p-value cutoff for a specific FDR
#' 
#' @param plierRes A PLIER result
#' @param fdr.cutoff The cross-validation significance cutoff for a pathway to be considered for naming
#' @export
getCutoff=function(plierRes,  fdr.cutoff=0.01){
max(plierRes$summary[plierRes$summary[,"FDR"]<=fdr.cutoff,"p-value"])
}


#' names latent variables according to their pathway useage (if any)
#' 
#' @param plierRes A PLIER result
#' @param top The number of pathway to use. Only the top pathway (one with the largest coefficient) is used by default
#' @param fdr.cutoff The cross-validation significance cutoff for a pathway to be considered for naming. If no pathways satisfy the cutoff the raw coefficients are used.
#' @param Choice of coef or AUC, whether LVs are named based on U coefficients or AUCs. Defualt: coef.
#' @export
nameB=function(plierRes, top=1, fdr.cutoff=0.01, use=c("coef", "AUC")){
use=match.arg(use, c("coef", "AUC"))
  names=vector("character",ncol(plierRes$U))
if(use=="coef"){
    Uuse=plierRes$U
}
  else{
    Uuse=plierRes$Uauc
  }
  if(!is.null(plierRes[["Up"]])){
    pval.cutoff=max(plierRes$summary[plierRes$summary[,5]<fdr.cutoff,4])
  
    Uuse[plierRes$Up>pval.cutoff]=0
  
  }
  else{
    warning("No p-values in PLIER object: using coefficients only")
  }
  mm=apply(Uuse,2,max)
  for(i in 1:ncol(plierRes$U)){
    if(mm[i]>0){
      names[i]=paste(i,names(sort(Uuse[,i],T))[1:top], sep=",")
    }
    else if(max(plierRes$U[,i])>0){
      names[i]=paste(i,names(sort(plierRes$U[,i],T))[1:top], sep=",")
    }
    else{
      names[i]=paste("LV",i)
    }
  }
  
  names
  
}

#' Computes the ridge pseudo-inverse of the prior information matrix. Used internally by PLIER but can be precomputed if running PLIER multiple times.
#' 
#' @param gsMat The prior information matrix. The genes have to match the gene expression data to be analyzed exactly (same genes and same order(
#' @param lambda The regularization paramter
#' @export
computeChat=function(gsMat, lambda=5){
  Chat=pinv.ridge(crossprod(gsMat,), lambda)%*%(t(gsMat))
}

#' @keywords internal
copyMat=function(mat, zero=F){
  matnew=matrix(nrow=nrow(mat), ncol=ncol(mat))
  rownames(matnew)=rownames(mat)
  colnames(matnew)=colnames(mat)
  if(zero)
    matnew[]=0
  matnew
}



#' @keywords internal
#' @param priorMat the real prior info matrix
#' @param priorMatcv the zeroed-out prior info matrix used for PLIER computations
#' 
crossVal=function(plierRes, data, priorMat, priorMatcv){
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=copyMat(plierRes$U,T)
  Up=copyMat(plierRes$U,T)
  Up[]=1
  for ( i in ii){
    
    iipath=which(plierRes$U[,i]>0)
    
    if (length(iipath) > 1){
    for(j in iipath){
      iiheldout=which((rowSums(priorMat[,iipath, drop=F])==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
      aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
      out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }}else{
      j <- iipath
      iiheldout=which((rowSums(matrix(priorMat[,iipath],ncol=1))==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
      aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
      out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
      }#else
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR") 
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}


#' @keywords internal
#' @param plierRes current PLIER result
#' @param data the input data 
#' @param priorMat the  prior info matrix

#' 
getAUC=function(plierRes, data, priorMat){
  Y=data
  B=plierRes$B
  Z=plierRes$Z
  Zcv=copyMat(Z)
  k=ncol(Z)
  L1=plierRes$L1
  L2=plierRes$L2
  for (i in 1:5){
    ii=(0:(floor(nrow(data)/5)-1))*5+i
    ii=ii[ii<=nrow(Z)]
    
    
    Bcv=solve(crossprod(Z[-ii,])+L2*diag(k))%*%t(Z[-ii,])%*%Y[-ii,]
    
    Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%solve(tcrossprod(Bcv)+L1*diag(k))
  }
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=copyMat(plierRes$U,T)
  Up=copyMat(plierRes$U,T)
  Up[]=1;
  for ( i in ii){
    
    iipath=which(plierRes$U[,i]>0)
    
    for(j in iipath){
      aucres=AUC(priorMat[,j], Zcv[,i])
      out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR") 
  
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}



#' Main PLIER function
#' 
#' @param data the data to be processed with genes in rows and samples in columns. Should be z-scored or set scale=T 
#' @param priorMat the binary prior information matrix with genes in rows and pathways/genesets in columns
#' @param svdres Pre-computed result of the svd decomposition for data
#' @param k The number of latent variables to return, leave as NULL to be set automatically using the num.pc "elbow" method
#' @param L1 L1 constant, leave as NULL to automatically select a value
#' @param L2 L2 constant, leave as NULL to automatically select a value
#' @param L3 L3 constant, leave as NULL to automatically select a value. Sparsity in U should be instead controlled by setting frac
#' @param frac The fraction of LVs that should have at least 1 prior inforamtion association, used to automatically set L3
#' @param max.iter Maximum number of iterations to perform
#' @param trace Display progress information
#' @param scale Z-score the data before processing
#' @param Chat A ridge inverse of priorMat, used to select active pathways, expensive to compute so can be precomputed when running PLIER multiple times
#' @param maxPath The maximum number of active pathways per latent variable
#' @param doCrossval Whether or not to do real cross-validation with held-out pathway genes. Alternatively, all gene annotations are used and only pseudo-crossvalidation is done. The latter option may be preferable if some pathways of interest have few genes. 
#' @param penalty.factor A vector equal to the number of columns in priorMat. Sets relative penalties for different pathway/geneset subsets. Lower penalties will make a pathway more likely to be used. Only the relative values matter. Internally rescaled.
#' @param glm_alpha Set the alpha for elastic-net
#' @param minGenes The minimum number of genes a pathway must have to be considered
#' @param tol Convergence threshold
#' @param seed Set the seed for pathway cross-validation
#' @param  allGenes Use all genes. By default only genes in the priorMat matrix are used.
#' @param rseed Set this option to use a random initialization, instead of SVD
#' @param pathwaySelection Pathways to be optimized with elstic-net penalty are preselected based on ridge regression results. "Complete" uses all top  pathways to fit individual LVs. "Fast" uses only the top pathways for the single LV in question.
#' @export

PLIER=function(data, priorMat,svdres=NULL, k=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=T, Chat=NULL, maxPath=10, doCrossval=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, pathwaySelection=c("complete", "fast")){
  
  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))

  
  
  
  
  if(scale){
    Y=rowNorm(data)
  }
  else{
    Y=data
  }
  
  if(nrow(priorMat)!=nrow(data) || !all(rownames(priorMat)==rownames(data))){
    if(!allGenes){
      cm=commonRows(data, priorMat)
      message(paste("Selecting common genes:", length(cm)))
      priorMat=priorMat[cm,]
      Y=Y[cm,]
    }
    else{
      extra.genes=setdiff(rownames(data), rownames(priorMat))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(priorMat))
      rownames(eMat)=extra.genes
      priorMat=rbind(priorMat, eMat)
      priorMat=priorMat[rownames(data),]
    }
    
  }
  numGenes=colSums(priorMat)
  
  heldOutGenes=list()
  iibad=which(numGenes<minGenes)
  priorMat[, iibad]=0
  message(paste("Removing", length(iibad), "pathways with too few genes"))
  if(doCrossval){
    
    
    priorMatCV=priorMat
    if(!is.null(seed))
      set.seed(seed)
    for(j in 1:ncol(priorMatCV)){
      
      iipos=which(priorMatCV[,j]>0)
      iiposs=sample(iipos, length(iipos)/5)
      priorMatCV[iiposs,j]=0
      heldOutGenes[[colnames(priorMat)[j]]]=rownames(priorMat)[iiposs]
      
    }
    C = priorMatCV
  }
  else{
    C=priorMat
  }
  
  nc=ncol(priorMat)
  ng=nrow(data)
  ns=ncol(data)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv.ridge(crossprod(C), 5)%*%(t(C))
  }
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y)){
    message("SVD V has the wrong number of columns")
    svdres=NULL
  }
  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }
    else{
      svdres=svd(Y)
    }
    message("Done")
  }
  if(is.null(k)){
    k=num.pc(svdres)*2
    k <- min(k, floor(ncol(Y)*0.9))
    message("k is set to ", k)
  }
  
  
  
  
  
  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }
  
  
  B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(svdres$d[1:k]))
  Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
  Z[Z<0]=0
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=apply(Z,2,sample)
  }
  
  
  
  
  
  
  U=matrix(0,nrow=ncol(C), ncol=k)
  
  
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  
  
  
  
  
  iter.full.start=iter.full=20
  
  curfrac=0
  nposlast=Inf
  npos=-Inf
  if(!is.null(L3)){
    L3.given=T
  }
  else{
    L3.given=F
  }
  
  for ( i in 1:max.iter){
    
    
    
    
    
    if(i>=iter.full.start){
      
      
      
      
      
      if(i==iter.full & !L3.given){ #update L3 to the target fraction
        Ulist=solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, target.frac = frac)
    U=Ulist$U
    L3=Ulist$L3
   message(paste("New L3 is", L3))
   iter.full=iter.full+iter.full.start
          }
      else{
      #HERE
      #solveU=function(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, target.frac=0.7, L3=NULL)
       
      U=solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, L3=L3)
      }
      curfrac=(npos<-sum(apply(U,2,max)>0))/k
      Z1=Y%*%t(B)
      Z2=L1*C%*%U
      ratio=median((Z2/Z1)[Z2>0&Z1>0])
      Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
    }
    
    else{
      Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    }
    
    
    Z[Z<0]=0
    
    
    oldB=B
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    
    
    
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    
    
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace & i >=iter.full.start){
      
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }
    else if (trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", i, " Bdiff is not decreasing"))
      break
    }
    
  }
  rownames(U)=colnames(priorMat)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  
  if(doCrossval){
    outAUC=crossVal(out, Y, priorMat, priorMatCV)
  }
  else{
    message("Not using cross-validation. AUCs and p-values may be over-optimistic")
    outAUC=getAUC(out, Y, priorMat)
  }
  out$withPrior=which(colSums(out$U)>0)
  out$Uauc=outAUC$Uauc
  out$Up=outAUC$Upval
  out$summary=outAUC$summary
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  
  rownames(out$B)=nameB(out)
  
  out
}


#' plot the U matrix from a PLIER decomposition
#' 
#' @param plierRes the result returned by PLIER
#' @param auc.cutoff the AUC cutoff for pathways to be displayed, increase to get a smaller subset of U
#' @param fdr.cutoff the significance cutoff for the pathway-LV association
#' @param indexCol restrict to a subset of the columns (LVs)
#' @param indexRow restrict to a subset of rows (pathways). Useful if only interested in pathways of a specific type
#' @param top the number of top pathways to discplay for each LV
#' @param sort.row do not custer the matrix but instead sort it to display the positive values close do the 'diagonal'
#' @param ... options to be passed to pheatmap
#' @export
plotU=function(plierRes, auc.cutoff=0.6, fdr.cutoff=0.05, indexCol=NULL, indexRow=NULL, top=3, sort.row=F,...){
  if(is.null(indexCol)){
    indexCol=1:ncol(plierRes$U)
  }
  if(is.null(indexRow)){
    indexRow=1:nrow(plierRes$U)
  }
  
  U=plierRes$U
  pval.cutoff=max(plierRes$summary[plierRes$summary[,5]<fdr.cutoff,4])
  U[plierRes$Uauc<auc.cutoff]=0
  U[plierRes$Up>pval.cutoff]=0
  
  U=U[indexRow, indexCol]
  for ( i in 1:ncol(U)){
    ct=sort(U[,i],T)[top]
    
    U[U[,i]<ct,i]=0
  }
  rownames(U)=strtrim(rownames(U), 30)
  if(sort.row){
    Utmp=sweep(sign(U), 2,1:ncol(U)*100,  "*")
    Um=apply(Utmp,1,max)
    show(Um[Um!=0])
    U=U[order(-Um),]
    plotMat(U, cluster_row=F,...)
  }
  else{
    plotMat(U, ...)
  }
}


#' visualize the top genes contributing to the LVs
#' 
#' @param plierRes the result returned by PLIER
#' @param data the data to be displayed in a heatmap, typically the z-scored input data (or some subset thereof)
#' @param priorMat the same gene by geneset binary matrix that was used to run PLIER
#' @param top the top number of genes to use
#' @param index the subset of LVs to display
#' @param regress remove the effect of all other LVs before plotting top genes, will take longer but can be useful to see distinct patterns in highly correlated LVs.
#' @param allLVs plot even the LVs that have no pathway association
#' @param ... Additional arguments to be passed to pheatmap, such as a column annotation data.frame (annotation_col). See ?pheatmap for details.
#' @export
plotTopZ=function(plierRes, data, priorMat, top=10, index=NULL, regress=F, allLVs=F,...){
  data=data[rownames(plierRes$Z),]
  priorMat=priorMat[rownames(plierRes$Z),]
  ii=which(colSums(plierRes$U)>0)
  if(!allLVs){
  if(! is.null(index)){
    ii=intersect(ii,index)
  }
  }
  else{
    ii=index
  }
  
  tmp=apply(-plierRes$Z[, ii, drop=F],2,rank)
  nn=character()
  nncol=character()
  nnpath=character()
  nnindex=double()
  for (i in 1:length(ii)){
    nn=c(nn,nntmp<-names(which(tmp[,i]<=top)))
    nncol=c(nncol, rep(rownames(plierRes$B)[ii[i]], length(nntmp)))
    nnpath=c(nnpath,rowSums(priorMat[nntmp,plierRes$U[,ii[i]]>0, drop=F])>0)
    nnindex=c(nnindex,rep(ii[i], length(nntmp)))
    
  }
  names(nncol)=nn
  nncol=strtrim(nncol, 30)
  
  nnrep=names(which(table(nn)>1))
  if(length(nnrep)>0){
    nnrep.im=match(nnrep,nn)
    nn=nn[-nnrep.im]
    nncol=nncol[-nnrep.im]
    nnpath=nnpath[-nnrep.im]
    nnindex=c(nnindex,rep(ii[i], length(nntmp)))
    
  }
  nnpath[nnpath=="TRUE"]="inPathway"
  nnpath[nnpath=="FALSE"]="notInPathway"
  
  nncol=as.data.frame(list(nncol,nnpath))
  
  names(nncol)=c("pathway", "present")
  ll=c(inPathway="black", notInPathway="beige")
  
  anncol=list(present=ll)
  toPlot=tscale(data[nn,])
  
  if(regress){
    for ( i in ii){
      gi=which(nnindex==i)
    #  toPlot[gi,]=toPlot[gi, ]-plierRes$Z[rownames(toPlot)[gi],-i ]%*%plierRes$B[-i,colnames(toPlot)]
      toPlot[gi, ] =resid(toPlot[gi,], model.matrix(~t(plierRes$B[-i, colnames(toPlot)])))                                            
    }
  }
  
  maxval=max(abs(toPlot))
  
  pheatmap(toPlot, breaks=seq(-maxval, maxval, length.out = 99),color=colorpanel(100, "green", "white", "red"),annotation_row=nncol, show_colnames = F, annotation_colors = anncol, ...)
}

#' @keywords  internal
colSumNorm=function(matrix, return.all=F){
  ss=sqrt(colSums(matrix^2))
  ss[ss<1e-16]=1
  if(!return.all){
    return(  sweep(matrix,2,ss, "/"))
  }
  else{
    return(list(mat=sweep(matrix,2,ss, "/"), ss=ss))
  }
}


#' returns the row names in common
#' @param data1 One matrix with gene rownames
#' @param data2 Another matrix with gene rownames
#' @export
commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}
#' z-score each row of the data
#' @param x gene-expression matrix, with genes in rows
#' @export
rowNorm=function(x){
  s=apply(x,1,sd)
  m=apply(x,1,mean);
  x=sweep(x,1,m)
  x=sweep(x,1,s,"/")
  x
}
#' estimates the number of 'significant' principle components for the SVD decomposition -- this is the minimum k for PLIER

#' @param  data the same data as to be used for PLIER (z-score recommended) or alternatively the result of an svd calculation 
#' @param method Either "eblow" (fast) or "permutation" (slower, but less heuristic)
#' @param B number of permutations
#' @param seed seed for reproducibility 
#' @export
num.pc = function (data, method="elbow", B = 20, seed = NULL) 
{

  method=match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  if((class(data)!="list")&(class(data)!="rsvd")){
    message("Computing svd")
  n <- ncol(data)
  m <- nrow(data)
  data=rowNorm(data)
  if(n<500){
    k=n
  }
  else{
    k=max(200,n/4)
  }
  if(k==n){
    uu <- svd(data)
  }
  else{
  set.seed(123456);uu <- rsvd(data,k, q=3)
  }
  }
  else if (!is.null(data[["d"]])){
    if(method=="permutation"){
      message("Original data is needed for permutation method.\nSetting method to elbow")
    method="elbow"
    }
  
    uu=data
    
  }
    
  
 
  if(method=="permutation"){
  #nn = min(c(n, m))
  dstat <- uu$d[1:k]^2/sum(uu$d[1:k]^2)
  dstat0 <- matrix(0, nrow = B, ncol = k)
  for (i in 1:B) {
    dat0 <- t(apply(data, 1, sample, replace = FALSE))
    if(k==n){
      uu0 <- svd(dat0)
    }
    else{
      set.seed(123456);
      uu0 <- rsvd(dat0,k, q=3)
    }
    dstat0[i, ] <- uu0$d[1:k]^2/sum(uu0$d[1:k]^2)
  }
  psv <- rep(1, k)
  for (i in 1:k) {
    psv[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:k) {
    psv[i] <- max(psv[(i - 1)], psv[i])
  }

  nsv <- sum(psv <= 0.1)
  }
  else if (method=="elbow"){
    x=smooth(xraw<-abs(diff(diff(uu$d))), twiceit = T)
 #plot(x)
 

 nsv=which(x<=quantile(x, 0.5))[1]+1
    
  }
  return(nsv)
}

#' @keywords internal
pinv.ridge=function (m, alpha = 0) 
{
  msvd = svd(m)
  if (length(msvd$d) == 0) {
    return(array(0, dim(m)[2:1]))
  }
  else {
    if (alpha > 0) {
      ss = (msvd$d^2) + alpha^2
      msvd$d = ss/msvd$d
    }
    out = msvd$v %*% (1/msvd$d * t(msvd$u))
    rownames(out) = rownames(m)
    colnames(out) = colnames(m)
    out
  }
}

#' generic function to plot the non-zero elements of a sparse matrix
#' @param  matrix the input matrix
#' @param scale rescale the columns to max=1
#' @param trim.names the maximum length of row and column lables
#' @param col.scale custom color scale
#' @param cutoff Sets values (both positive and negative) bellow this number to 0
#' @param ... additional argumend to be passed to pheatmap
#' @export
plotMat=function(matrix,  scale=T, trim.names=50, cutoff=NULL,col.scale=NULL,...){
  
  if(! is.null(trim.names)){
    rownames(matrix)=strtrim(rownames(matrix), trim.names)
    colnames(matrix)=strtrim(colnames(matrix), trim.names)
  }
  if(!is.null(cutoff)){
    matrix[abs(matrix)<cutoff]=0
  }
  matrix=matrix[iirow<-rowSums(abs(matrix))>0,]
  matrix=matrix[,iicol<-colSums(abs(matrix))>0]
  
  mydistBin=function(x){dist(abs(sign(x)))}

  mydistBin=function(x){dist(abs(sign(x)))}
  
  if(scale){
    
    aa=apply(abs(matrix),2, max)
    aa[aa==0]=1
    
    matrix=sweep(matrix,2,aa, "/")
    
  }
  if (min(matrix)<0)
    mycol= c("grey90",colorRampPalette(rev(brewer.pal(n = 7, "RdYlBu")))(100))
  else
    mycol=c("white",colorRampPalette(rev(brewer.pal(n = 11, name =  "PRGn"))[7:11])(100))
  if(!is.null(col.scale)){
    mycol=colscale
  }
 
  pheatmap(matrix,color =mycol , clustering_callback = function(h,d){hclust(mydistBin(d), method = "single")}, ...)

  return(invisible(list(iirow=iirow, iicol=iicol)))
}

#' @keywords internal
tscale=function(mat){
  t(scale(t(mat)))
}
#' visualize the top genes contributing to the LVs similarily to \code{\link{plotTopZ}}. However in this case all the pathways contributing to each LV are show seperatly. Useful for seeing pathway usage for a single LV or understading the differences between two closely related LVs
#' 
#' @param plierRes the result returned by PLIER
#' @param data the data to be displayed in a heatmap, typically the z-scored input data (or some subset thereof)
#' @param priorMat the same gene by geneset binary matrix that was used to run PLIER
#' @param top the top number of genes to use
#' @param index the subset of LVs to display
#' @param regress remove the effect of all other LVs before plotting top genes, will take longer but can be useful to see distinct patterns in highly correlated genes.
#' @param fdr.cutoff Significance cutoff for a pathway to be plotted
#' @param ... Additional arguments to be passed to pheatmap, such as a column annotation data.frame
#' @export
#' 
#' 
plotTopZallPath=function (plierRes, data, priorMat, top = 10, index = NULL, regress = F, 
                          fdr.cutoff = 0.2, ...) 
{
  pval.cutoff = max(plierRes$summary[plierRes$summary[, 5] < 
                                       fdr.cutoff, 4])
  ii = which(colSums(plierRes$U) > 0)
  if (!is.null(index)) {
    ii = intersect(ii, index)
  }
  tmp = apply(-plierRes$Z[, ii, drop = F], 2, rank)
  nn = character()
  nncol = character()
  nnpath = character()
  nnindex = double()
  Ustrict = plierRes$U
  Ustrict[plierRes$Up > pval.cutoff] = 0
  pathsUsed = which(rowSums(Ustrict[, index, drop = F]) > 0)
  pathMat = matrix(nrow = 0, ncol = length(pathsUsed))

    colnames(pathMat) = strtrim(names(pathsUsed), 40)
  
  
  #  colnames(pathMat) = wrapString(names(pathsUsed), 40)

  for (i in 1:length(ii)) {
    nn = c(nn, nntmp <- names(which(tmp[, i] <= top)))
    nncol = c(nncol, rep(rownames(plierRes$U)[which(thispath <- plierRes$U[, 
                                                                           ii[i]] == max(plierRes$U[, ii[i]]))], length(nntmp)))
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
    pathMat = rbind(pathMat, priorMat[nntmp, pathsUsed, drop=F])
  }
  
  if(sum(colSums(pathMat)>1)>0){
    pathMat = pathMat[, colSums(pathMat) > 0]
    pathsUsed=colnames(pathMat)
  }
  pathMat = as.data.frame(pathMat)

  
  pathMat = apply(pathMat, 2, as.factor)
  
  names(nncol) = nn
  nncol = strtrim(nncol, 40)
  nnrep = names(which(table(nn) > 1))
  ll = list(inPathway = "black", notInPathway = "beige")
  ll2 = list()
  for (i in 1:length(pathsUsed)) {
    ll2[[i]] = c("black", "beige")
    names(ll2[[i]]) = c("1", "0")
  }
  names(ll2) = colnames(pathMat)
  
  anncol = ll2
 
  rr = max(range(tscale(data[nn, ])))
  bb = seq(-rr, rr, length.out = 100)
  toPlot = data[nn, ]
  if (regress) {
    for (i in ii) {
      gi = which(nnindex == i)
 
    #  toPlot[gi, ] = toPlot[gi, ] - plierRes$Z[rownames(toPlot)[gi], 
#                                               -i] %*% plierRes$B[-i, colnames(toPlot)]
#                                               
      toPlot[gi, ] =resid(toPlot[gi,], model.matrix(~t(plierRes$B[-i, colnames(toPlot)])))                                            
    }
  }
  pheatmap(tscale(toPlot), breaks = bb, color = colorpanel(101, 
                                                           "green", "white", "red"), annotation_row = as.data.frame(pathMat
                                                           ), annotation_legend = F, show_colnames = F, annotation_colors = anncol, 
           clustering_callback = function(h, d) {
             hclust(dist(d), method = "complete")
           }, ...)
}
#' @keywords internal
nonEstimable=function (x) 
{
  x = as.matrix(x)
  p = ncol(x)
  QR = qr(x)
  if (QR$rank < p) {
    n = colnames(x)
    if (is.null(n)) 
      n = as.character(1:p)
    notest = n[QR$pivot[(QR$rank + 1):p]]
    blank = notest == ""
    if (any(blank)) 
      notest[blank] = as.character(((QR$rank + 1):p)[blank])
    return(notest)
  }
  else {
    return(NULL)
  }
}
#' @keywords internal
resid=function(dat, lab, useMean=T){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  ne = nonEstimable(mod)
  if (!is.null(ne)){ 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
    mod=mod[, -match(ne, colnames(mod))]
  }
  
  n=dim(dat)[2]
  Id=diag(n)
  out=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                 t(mod))
  colnames(out)=colnames(dat)
  if (useMean){
    out=sweep(out,1,apply(dat,1,mean), "+")  
  }
  
  out
}




#' @keywords internal
scad=function(x, lambda,a=3.7){
  
  iip=which(abs(x)>2*lambda & abs(x)<a*lambda)
  
  iin=which(abs(x)<=2*lambda)
  
  x[iin]=pmax(0, abs(x[iin])-lambda)*sign(x[iin])
  x[iip]=((a-1)*x[iip]-sign(x[iip])*a*lambda)/(a-2)
  
  x
}
#' @keywords internal
quicksoft=function (x, d) {
 return(sign(x) * pmax(0, abs(x) - d))
}


#' @keywords internal
scadZ=function(Z, ii=1:ncol(Z), lambda){
  Zn=colSumNorm(Z, return.all = T)
  Zt=Z
  
  #Zt[,ii]=apply(Zn$mat[,ii], 2, function(x){scad(x,lambda)})
  #Zt[,ii]=sweep(Zt[, ii],2, Zn$ss[ii], "*")
  Zt[,ii]=apply(Z[,ii], 2, function(x){scad(x, lambda)})
  Zt
}

#' @keywords internal
softZ=function(Z, ii=1:ncol(Z), lambda){
  Zn=colSumNorm(Z, return.all = T)
  Zt=Z
  
  #Zt[,ii]=apply(Zn$mat[,ii], 2, function(x){quicksoft(x,lambda)})
  #Zt[,ii]=sweep(Zt[, ii],2, Zn$ss[ii], "*")
   Zt[,ii]=apply(Z[,ii], 2, function(x){quicksoft(x, lambda)})
  Zt
}
#' @keywords internal
getEnrichmentVals=function(plierRes, pathwayMat, ngenes=50,auc.cutoff=0.7, fdr.cutoff=0.01){
  pathwayMat=pathwayMat[rownames(plierRes$Z), rownames(plierRes$U)]
  Uuse=plierRes$U
  Uuse[plierRes$Uauc<auc.cutoff]=0
  Uuse[plierRes$Up>getCutoff(plierRes, fdr.cutoff)]=0
  inpath=intop=double(ncol(plierRes$Z))
  
  for(i in 1:ncol(plierRes$Z)){
    iipath=which(Uuse[,i]>0)
    if(length(iipath)>0){
      pathGenes=names(which(rowSums(pathwayMat[,iipath, drop=F])>0))
      topGenes=names(sort(plierRes$Z[,i], T)[1:ngenes])
      pathGenesInt=intersect(pathGenes, topGenes)
      inpath[i]=length(pathGenes)
      intop[i]=length(pathGenesInt)
    }}
  return(cbind(intop/inpath,intop, inpath))
}



#' sparse PLIER function (experimental)
#' 
#' @param data the data to be processed with genes in rows and samples in columns. Should be z-scored or set scale=T 
#' @param priorMat the binary prior information matrix with genes in rows and pathways/genesets in columns
#' @param svdres Pre-computed result of the svd decomposition for data
#' @param k The number of latent variables to return, leave as NULL to be set automatically using the num.pc "elbow" method
#' @param L1 L1 constant, leave as NULL to automatically select a value
#' @param L2 L2 constant, leave as NULL to automatically select a value
#' @param L3 L3 constant, leave as NULL to automatically select a value. Sparsity in U should be instead controlled by setting frac
#' @param frac The fraction of LVs that should have at least 1 prior inforamtion association, used to automatically set L3
#' @param max.iter Maximum number of iterations to perform
#' @param trace Display progress information
#' @param scale Z-score the data before processing
#' @param Chat A ridge inverse of priorMat, used to select active pathways, expensive to compute so can be precomputed when running PLIER multiple times
#' @param maxPath The maximum number of active pathways per latent variable
#' @param doCrossval Whether or not to do real cross-validation with held-out pathway genes. Alternatively, all gene annotations are used and only pseudo-crossvalidation is done. The latter option may be preferable if some pathways of interest have few genes. 
#' @param penalty.factor A vector equal to the number of columns in priorMat. Sets relative penalties for different pathway/geneset subsets. Lower penalties will make a pathway more likely to be used. Only the relative values matter. Internally rescaled.
#' @param glm_alpha Set the alpha for elastic-net
#' @param minGenes The minimum number of genes a pathway must have to be considered
#' @param tol Convergence threshold
#' @param seed Set the seed for pathway cross-validation
#' @param  allGenes Use all genes. By default only genes in the priorMat matrix are used.
#' @param rseed Set this option to use a random initialization, instead of SVD
#' @param pathwaySelection Pathways to be optimized with elstic-net penalty are preselected based on ridge regression results. "Complete" uses all top  pathways to fit individual LVs. "Fast" uses only the top pathways for the single LV in question.
#' @param sparseL the lambda for sparsity on Z, default 0.02
#' @param sparseType sparsity inducing penalty to use (SCAD or L1)
#' @export
PLIERsparse=function(data, priorMat,svdres=NULL, k=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=T, Chat=NULL, maxPath=10, doCrossval=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, pathwaySelection=c("complete", "fast"), sparseL=0.01, sparseType="SCAD"){
  
  sparseType=match.arg(sparseType, c("SCAD", "L1"))
  
  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))
  #Ur is the ranked matrix of pathway relevance
  solveU=function(Z, Ur,C,  L3, penalty.factor, glm_alpha){
    
    ii=which(apply(Ur,1,min)<=maxPath)
    
    U=copyMat(Ur)
    U[]=0
    
    for (j in 1:ncol(U)){
      if(pathwaySelection=="fast"){
        selection=which(Ur[,j]<=maxPath)
      }
      else{
        selection=ii
      }
      #  if(conditional){
      #    iigenes=which(Z[,j]>0)
      #  }
      #  else{
      iigenes=1:nrow(Z)
      #  }
      Zr=rank(-Z[iigenes])
    
      tmp=glmnet(y=Z[iigenes,j],  x=priorMat[iigenes,selection], alpha=glm_alpha, lambda=L3, lower.limits = 0, penalty.factor = penalty.factor[selection])
      U[selection,j]=as.numeric(tmp$beta)
    }
    
    return(U)
  }
  
  
  
  
  if(scale){
    Y=rowNorm(data)
  }
  else{
    Y=data
  }
  
  if(nrow(priorMat)!=nrow(data) || !all(rownames(priorMat)==rownames(data))){
    if(!allGenes){
      cm=commonRows(data, priorMat)
      message(paste("Selecting common genes:", length(cm)))
      priorMat=priorMat[cm,]
      Y=Y[cm,]
    }
    else{
      extra.genes=setdiff(rownames(data), rownames(priorMat))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(priorMat))
      rownames(eMat)=extra.genes
      priorMat=rbind(priorMat, eMat)
      priorMat=priorMat[rownames(data),]
    }
    
  }
  numGenes=colSums(priorMat)
  
  heldOutGenes=list()
  iibad=which(numGenes<minGenes)
  priorMat[, iibad]=0
  message(paste("Removing", length(iibad), "pathways with too few genes"))
  if(doCrossval){
    
    
    priorMatCV=priorMat
    if(!is.null(seed))
      set.seed(seed)
    for(j in 1:ncol(priorMatCV)){
      iipos=which(priorMatCV[,j]>0)
      iiposs=sample(iipos, length(iipos)/5)
      priorMatCV[iiposs,j]=0
      heldOutGenes[[colnames(priorMat)[j]]]=rownames(priorMat)[iiposs]
      
    }
    C = priorMatCV
  }
  else{
    C=priorMat
  }
  
  nc=ncol(priorMat)
  ng=nrow(data)
  ns=ncol(data)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat) || (pathwaySelection=="fast"&&doCrossval)){
    Cp=crossprod(C)
    Chat=pinv.ridge(crossprod(C), 5)%*%(t(C))
  }
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  
  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }
    else{
      svdres=svd(Y)
    }
    message("Done")
  }
  if(is.null(k)){
    k=num.pc(svdres)*2
    message("k is set to ", k)
  }
  if(nrow(svdres$u)!=nrow(Y)){
    message("SVD U has the wrong number of rows")
    
    if(!is.null(rownames(svdres$u))){
      message("Selecting via rownames of U")
      Z=svdres$u[rownames(Y),1:k]
    }
    else{
      message("No rownames for svdres$u: Recomputing SVD")
      if(ns>500){
        if(!is.null(seed))
          set.seed(seed)
        set.seed(123456);svdres=rsvd(Y, k, q=3)
      }
      else{
        svdres=svd(Y)
      }
      
      message("Done")
    }
  }
  else{
    Z=svdres$u[, 1:k]
  }
  
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=t(apply(Z,2,sample))
  }
  
  
  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }
  
  
  
  B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(svdres$d[1:k]))
  U=matrix(0,nrow=ncol(C), ncol=k)
  
  show(dim(B))
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  
  Z[Z<0]=0
  
  
  
  iter.full.start=iter.full=20
  
  curfrac=0
  nposlast=Inf
  npos=-Inf
  if(!is.null(L3)){
    L3.given=T
  }
  else{
    L3.given=F
  }
  
  for ( i in 1:max.iter){
    
    
    
    
    
    if(i>=iter.full.start){
      
      #Compute Us
      Us=Chat%*%colSumNorm(Z)
      
      Us[Us<0]=0
      Us=apply(-Us,2,rank)
      
      #    ii=which(apply(Us,1,min)<=maxPath)
      
      if(i==iter.full & !L3.given){
        
        
        message(paste0("Updating L3, current fraction= ", round(curfrac,4), ", target=", frac))
        biter=0
        
        
        if(abs(frac-curfrac)>1/k){
          #set up the limits
          if(curfrac>frac){
            #increase penatly
            if(is.null(L3)){
              L3_1=0.000001
              L3_2=1
            }
            else{
              L3_1=L3
              L3_2=L3*100
            }
          }
          else{
            #decrease
            if(is.null(L3)){
              
              L3_1=0.000001
              L3_2=1
            }
            else{
              L3_1=L3/100
              L3_2=L3
            }
          }
          
          
          while (biter < 150&(biter<1|abs(frac-curfrac)>1/k|npos==0)){
            
            U=solveU(Z, Us, C,  L3=(L3use<-(L3_1+L3_2)/2), penalty.factor, glm_alpha)
            
            nposlast=npos
            curfrac=(npos<-sum(apply(U,2,max)>0))/k
            if(T){
              message(paste0(npos, " positive columns at L3=", round(L3use,6)))
            }
            if(curfrac>frac){  #increase penatly
              #check if the limits have been reached
              if((L3_2-L3_1)<1e-7){
                L3_2=L3_2*100 
              }
              L3_1=(L3_1+L3_2)/2
            }
            else{#decrease
              if((L3_2-L3_1)<1e-7){
                L3_1=L3_1/100 
              }
              L3_2=(L3_1+L3_2)/2
              
            }
            
            biter=biter+1
            #show(c(npos, nposlast, frac, curfrac, abs(frac-curfrac), 1/k))
          }
          L3=L3use
          if(trace){
            message(paste0("L3 is set to ", round(L3, 6), " in ", biter, " iterations"))
          }
        }
        else{
          message("L3 not changed")
        }
        iter.full=iter.full+iter.full.start
      }
      
      #find the active pathways
      #  ii=which(apply(Us,1,min)<=maxPath)
      U=solveU(Z, Us, C, L3, penalty.factor, glm_alpha)
      curfrac=(npos<-sum(apply(U,2,max)>0))/k
      Z1=Y%*%t(B)
      Z2=L1*C%*%U
      Z2[Z==0]=0
      ratio=median((Z2/Z1)[Z2>0&Z1>0])
      Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
    }
    
    else{
      Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    }
    
    Zraw=Z
    Z[Z<0]=0
    iipath=which(colSums(U)>0)
    
    if(length(iipath)>0){
      if(sparseType=="SCAD"){
        Z=scadZ(Z, iipath, lambda = sparseL)
      }
      else{
        Z=softZ(Z, iipath, lambda = sparseL)
      }
    }
    
    
    oldB=B
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    
    
    
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    
    
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace & i >=iter.full.start){
      
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }
    else if (trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", i, " Bdiff is not decreasing"))
      break
    }
    
  }
  rownames(U)=colnames(priorMat)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C,  L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes, sparseL=sparseL, sparseType=sparseType)
  
  if(doCrossval){
    outAUC=crossVal(out, Y, priorMat, priorMatCV)
  }
  else{
    message("Not using cross-validation. AUCs and p-values may be over-optimistic")
    outAUC=getAUC(out, Y, priorMat)
  }
  out$withPrior=which(colSums(out$U)>0)
  out$Uauc=outAUC$Uauc
  out$Up=outAUC$Upval
  out$summary=outAUC$summary
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  out$Zraw=Zraw
  rownames(out$B)=nameB(out)
  out
}


rotateSVD=function(svdres){
  sumposu=apply(svdres$u>0,2, sum)
  sumposv=apply(svdres$v>0,2, sum)
  
  sumpos=sumposu+sumposv
  iibad=which(sumpos==0)
  if(length(iibad)>0){
    message("All negative components found")
    stop()
  }
  for(i in 1:ncol(svdres$u)){
    if(sumposu[i]==0 | sumposv[i]==0){
      message("flipping")
      svdres$u[,i]=-svdres$u[,i] 
      svdres$v[,i]=-svdres$v[,i] 
    }
  }
  svdres
}#rotateSVD

simpleDecomp=function(Y, k,svdres=NULL, L1=NULL, L2=NULL,
         Zpos=T,max.iter=200, tol=5e-6, trace=F,
         rseed=NULL, B=NULL, scale=1, pos.adj=3, cutoff=0){
  
  
  
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=Inf
  BdiffTrace=double()
  BdiffCount=0
  message("****")
  
  if(is.null(svdres)){
    
    message("Computing SVD")
    set.seed(123);svdres=rsvd(Y, k = k) 
    
    svdres=rotateSVD(svdres)
    
    #  show(svdres$d[k])
  }
  
  
  if(is.null(L1)){
    L1=svdres$d[k]*scale
    if(!is.null(pos.adj)){
      L1=L1/pos.adj
    }
    
  }
  if(is.null(L2)){
    L2=svdres$d[k]*scale
  }
  #    L1=svdres$d[k]/2*scale
  print(paste0("L1 is set to ",L1))
  print(paste0("L2 is set to ",L2))
  
  if(is.null(B)){
    #initialize B with svd
    message("Init")
    B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(sqrt(svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k])
  }
  else{
    message("B given")
  }
  
  
  
  
  
  if (!is.null(rseed)) {
    message("using random start")
    set.seed(rseed)
    B = t(apply(B, 1, sample))
    
  }
  
  
  round2=function(x){signif(x,4)}
  
  
  
  
  for ( i in 1:max.iter){
    #main loop    
    Zraw=Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    
    if(Zpos){
      Z[Z<cutoff]=0
    }
    
    oldB=B
    
    
    B=solve(t(Z)%*%Z+L2*diag(k))%*%(t(Z)%*%Y)
    
    
    #update error
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z)^2)*L1+sum(B^2)*L2
    if(trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }
    
  }
  rownames(B)=colnames(Z)=paste("LV",1:k)
  return(list(B=B, Z=Z, Zraw=Zraw, L1=L1, L2=L2))
}#simpleDecomp
