library(RColorBrewer)
library(gplots)
library(pheatmap)
library(glmnet)

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


AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  if(posn>0&negn>0){
    res=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE);
    myres=list()
    myres$low=res$conf.int[1]
    myres$high=res$conf.int[2]
    myres$auc=(res$statistic)/(posn*negn)
    myres$pval=res$p.value
  }
  else{
    myres$auc=0.5
    mures$pval=NA
  }
  return(myres)
}


nameB=function(pres, top=1){
  names=vector("character",ncol(pres$U))
  mm=apply(pres$U,2,max)
  for(i in 1:ncol(pres$U)){
    if(mm[i]>0){
      names[i]=paste(i,names(sort(pres$U[,i],T))[1:top], sep=",")
    }
    else{
      names[i]=paste("LV",i)
    }
  }
  
  names
  
}


computeChat=function(gsMat, lambda=5){
  Chat=pinv.ridge(crossprod(gsMat,), lambda)%*%(t(gsMat))
}


copyMat=function(mat, zero=F){
  matnew=matrix(nrow=nrow(mat), ncol=ncol(mat))
  rownames(matnew)=rownames(mat)
  colnames(matnew)=colnames(mat)
  if(zero)
    matnew[]=0
  matnew
}


crossVal=function(plierRes, data, priorMat){
  Y=data
  B=plierRes$B
  Z=plierRes$Z
  Zcv=copyMat(Z)
  
  for (i in 1:5){
    ii=(0:(floor(nrow(data)/5)-1))*5+i
    ii=ii[ii<=nrow(Z)]
    
    Bcv=pinv.ridge(crossprod(Z[-ii,]))%*%t(Z[-ii,])%*%Y[-ii,]
    
    Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%pinv.ridge(tcrossprod(Bcv))
  }
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=copyMat(plierRes$U,T)
  Up=copyMat(plierRes$U,T)
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
  colnames(out)=c("pathway", "LV index", "AUC", "p-value") 
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}


PLIER=function(data, priorMat,svdres=NULL, k, L1=NULL, L2=NULL, L3=NULL,  frac=0.5, posU=T, posZ=T, max.iter=200, trace=F, scale=F, Chat=NULL, maxPath=10, computeAUC=T){
  
  glm_alpha=0.9
  if(scale){
    Y=rowNorm(data)
  }
  else{
    Y=data
  }
  
  C=priorMat
  
  
  nc=ncol(priorMat)
  ns=nrow(data)
  ng=ncol(data)
  U=matrix(0,nrow=nc, ncol=k )
  Z=matrix(0, nrow=ng, ncol=k)
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv.ridge(Cp, 5)%*%(t(C))
  }
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  if(is.null(svdres)){
    message("Computing SVD")
    svdres=svd(Y)
    message("Done")
  }
  if(nrow(svdres$u)!=nrow(data)){
    warning("SVD U has the wrong number of rows")
    
    if(!is.null(rownames(svdres$u))){
      warning("Selecting via rownames")
      Z=svdres$u[rownames(Y),1:k]
    }
    else{
      message("Computing SVD")
      svdres=svd(Y)
      message("Done")
    }
  }else{
    Z=svdres$u[, 1:k]
  }
  
  if(is.null(L2)){
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }
  
  
  
  B=t(svdres$v[, 1:k]%*%diag(svdres$d[1:k]))
  
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  Z[Z<0]=0
  for ( i in 1:max.iter){
    
    
    
    Us=Chat%*%colSumNorm(Z)
    
    
    if(i==1 & is.null(L3)){
      
      L3=quantile(apply(Us,2,max), 1-frac)
      print(paste0("L3 is set to ",L3))
    }
    
    Us[Us<0]=0
    Us=apply(-Us,2,rank)
    ii=which(apply(Us,1,min)<=maxPath)
    
    if(trace)
      message(paste("there are ", length(ii), "active pathways"))
    U[]=0
    
    for (j in 1:ncol(U)){
      tmp=glmnet(y=Z[,j], x=priorMat[, ii], alpha=glm_alpha, lambda=L3, lower.limits = 0)
      U[ii,j]=as.numeric(tmp$beta)
      
    }
    
    
    Z1=Y%*%t(B)
    Z2=L1*C%*%U
    
    ratio=median((Z2/Z1)[Z2>0&Z1>0])
    
    
    Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
    #  Z=colSumNorm(Z)
    
    if(posZ){
      Z[Z<0]=0
    } 
    
    oldB=B
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    
    
    
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    
    
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace){
      
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }
    
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<5e-6){
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
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, numActPath=length(ii), L1=L1, L2=L2, L3=L3)
  if(computeAUC){
    crossval=crossVal(out, Y,C)
    rownames(out$B)=nameB(out)
    out$withPrior=which(colSums(out$U)>0)
    out$Uauc=crossval$Uauc
    out$Up=crossval$Upval
    out$summary=crossval$summary
    tt=apply(out$Uauc,2,max)
    message(paste("There are", sum(tt>0.75), " LVs with AUC>0.75"))
  }
  out
}


plotU=function(plierRes, auc.cutoff=0.6, pval.cutoff=1e-05, index){
  U=plierRes$U
  U[plierRes$Uauc<auc.cutoff]=0
  U[plierRes$Up>pval.cutoff]=0
  plotMat(U[, index])
  
}


plotTopZ=function(resP, data, priorInfo, top=10, index=NULL){
  
  ii=which(colSums(resP$U)>0)
  if(! is.null(index)){
    ii=intersect(ii,index)
  }
  tmp=apply(-resP$Z[, ii, drop=F],2,rank)
  nn=character()
  nncol=character()
  nnpath=character()
  for (i in 1:length(ii)){
    nn=c(nn,nntmp<-names(which(tmp[,i]<=top)))
    nncol=c(nncol, rep(rownames(resP$U)[which(thispath<-resP$U[,ii[i]]==max(resP$U[,ii[i]]))], length(nntmp)))
    nnpath=c(nnpath,rowSums(priorInfo[nntmp,resP$U[,ii[i]]>0, drop=F])>0)
  }
  names(nncol)=nn
  nncol=strtrim(nncol, 30)
  
  nnrep=names(which(table(nn)>1))
  if(length(nnrep)>0){
    nnrep.im=match(nnrep,nn)
    nn=nn[-nnrep.im]
    nncol=nncol[-nnrep.im]
    nnpath=nnpath[-nnrep.im]
  }
  nnpath[nnpath=="TRUE"]="inPathway"
  nnpath[nnpath=="FALSE"]="notInPathway"
  
  nncol=as.data.frame(list(nncol,nnpath))
  
  names(nncol)=c("pathway", "present")
  ll=c(inPathway="black", notInPathway="beige")
  
  anncol=list(present=ll)
  
  pheatmap(tscale(data[nn,]), color=colorpanel(100, "green", "white", "red"),annotation_row=nncol, show_colnames = F, annotation_colors = anncol)
}


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


commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}

rowNorm=function(x){
  s=apply(x,1,sd)
  m=apply(x,1,mean);
  x=sweep(x,1,m)
  x=sweep(x,1,s,"/")
  x
}

num.pc = function (dat, B = 20, seed = NULL) 
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  n <- ncol(dat)
  m <- nrow(dat)
  uu <- svd(dat)
  nn = min(c(n, m))
  dstat <- uu$d[1:nn]^2/sum(uu$d[1:nn]^2)
  dstat0 <- matrix(0, nrow = B, ncol = nn)
  for (i in 1:B) {
    dat0 <- t(apply(dat, 1, sample, replace = FALSE))
    uu0 <- svd(dat0)
    dstat0[i, ] <- uu0$d[1:nn]^2/sum(uu0$d[1:nn]^2)
  }
  psv <- rep(1, nn)
  for (i in 1:nn) {
    psv[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:nn) {
    psv[i] <- max(psv[(i - 1)], psv[i])
  }
  show(psv)
  nsv <- sum(psv <= 0.1)
  return(as.numeric(list(n.sv = nsv)))
}


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


plotMat=function(matrix,  scale=T, trim.names=30,...){
  if(! is.null(trim.names)){
    rownames(matrix)=strtrim(rownames(matrix), trim.names)
  }
  matrix=matrix[iirow<-rowSums(abs(matrix))>0,]
  matrix=matrix[,iicol<-colSums(abs(matrix))>0]
  mydist=function(x){dist(abs(sign(x)))}
  
  if(scale){
    aa=apply(abs(matrix),2, max)
    aa[aa==0]=1  
    matrix=sweep(matrix,2,aa, "/")
    
  }
  
  pheatmap(matrix,color = c("white",colorRampPalette((brewer.pal(n = 7, name =  "BuPu")[2:7]))(100)), clustering_callback = function(h,d){hclust(mydist(d))}, ...)
  return(invisible(list(iirow=iirow, iicol=iicol)))
}


tscale=function(mat){
  t(scale(t(mat)))
}
