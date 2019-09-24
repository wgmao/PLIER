library(MASS)
library(glmnet)


tscale=function(mat){
  t(scale(t(mat)))
}#tscale

normF <- function(x){
  return(sum(x^2))
}#normF


rowNorm=function(x){
  s=apply(x,1,sd)
  m=apply(x,1,mean);
  x=sweep(x,1,m)
  x=sweep(x,1,s,"/")
  x
}#rowNorm

num.pc = function (data, method="elbow", B = 20, seed = NULL){
  
  method=match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  if(class(data)!="list"){
    message("Computing svd")
    n <- ncol(data)
    m <- nrow(data)
    data=rowNorm(data)
    if(n<500){
      k=n
    }else{
      k=max(200,n/4)
    }
    if(k==n){
      uu <- svd(data)
    }else{
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
    nn = min(c(n, m))
    dstat <- uu$d[1:nn]^2/sum(uu$d[1:nn]^2)
    dstat0 <- matrix(0, nrow = B, ncol = nn)
    for (i in 1:B) {
      dat0 <- t(apply(data, 1, sample, replace = FALSE))
      if(k==n){
        uu0 <- svd(dat0)
      }else{
        set.seed(123456);uu0 <- rsvd(dat0,k, q=3)
      }
      dstat0[i, ] <- uu0$d[1:nn]^2/sum(uu0$d[1:nn]^2)
    }
    psv <- rep(1, nn)
    for (i in 1:nn) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:nn) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }
    
    nsv <- sum(psv <= 0.1)
  }else if (method=="elbow"){
    x=smooth(xraw<-abs(diff(diff(uu$d))), twiceit = T)
    #plot(x)    
    nsv=which(x<quantile(x, 0.5))[1]+1
  }
  return(nsv)
}#num.pc


solveU=function(Z,  Chat, priorMat, penalty.factor, glm_alpha, target.frac=0.7){
  lambdas=exp(seq(-12,-4,0.125))
  
  Ur=Chat%*%Z #get U by OLS
  Ur=apply(-Ur,2,rank) #rank
  
  Urm=apply(Ur,1,min)
  ii=1:length(Z)
  results=list()
  lMat=matrix(nrow=length(lambdas), ncol=ncol(Z))
  for(i in 1:ncol(Z)){
    if(pathwaySelection=="fast"){
      iip=which(Ur[,i]<=maxPath)
    }else{
      iip=which(Urm<=maxPath)
    }#else
    gres=glmnet(y=Z[,i], x=priorMat[,iip], penalty.factor = penalty.factor, alpha=glm_alpha, lower.limits=0, lambda = lambdas,
                intercept=T,  standardize=T, )
    #plot(gres)
    gres$iip=iip
    lMat[,i]=colSums(gres$beta>0)
    results[[i]]=gres
  }
  fracs=rowMeans(lMat>0)
  iibest=which.min(abs(target.frac-fracs))
  iibest
  
  U=matrix(0,nrow=ncol(priorMat), ncol=ncol(Z))
  for(i in 1:ncol(Z)){
    U[results[[i]]$iip,i]=results[[i]]$beta[,iibest]
  }#for i
  rownames(U)=colnames(priorMat)
  colnames(U)=1:ncol(Z)
  U  
}#solveU




PLIERgrp=function(data, priorMat,svdres=NULL, k=NULL, grp=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=T, Chat=NULL, maxPath=10, doCrossval=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, pathwaySelection=c( "fast","complete")){
  
  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))
  #Ur is the ranked matrix of pathway relevance
  
  if(scale){
    Y=rowNorm(data)
  }else{
    Y=data
  }#if
  
  if(nrow(priorMat)!=nrow(data) || !all(rownames(priorMat)==rownames(data))){
    if(!allGenes){
      cm=commonRows(data, priorMat)
      message(paste("Selecting common genes:", length(cm)))
      priorMat=priorMat[cm,]
      Y=Y[cm,]
    }else{
      extra.genes=setdiff(rownames(data), rownames(priorMat))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(priorMat))
      rownames(eMat)=extra.genes
      priorMat=rbind(priorMat, eMat)
      priorMat=priorMat[rownames(data),]
    }#else
  }#if(nrow(priorMat)!=nrow(data) || !all(rownames(priorMat)==rownames(data)))
  
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
    }#for j
    C = priorMatCV
  }# if doCrossval
  else{
    C=priorMat
  }#else
  
  nc=ncol(priorMat)
  ng=nrow(data)
  ns=ncol(data)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv.ridge(crossprod(C), 5)%*%(t(C))
  }#if(is.null(Chat))
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y)){
    message("SVD V has the wrong number of columns")
    svdres=NULL
  }#if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y))
  
  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }else{
      svdres=svd(Y)
    }#else
    message("Done")
  }#if(is.null(svdres))
  
  if(is.null(k)){
    k=num.pc(svdres)*2
    message("k is set to ", k)
  }#if(is.null(k))
  
  
  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }#if(is.null(L2))
  
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }#if(is.null(L1))
  
  
  B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(svdres$d[1:k]))
  Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
  
  Zgrp_index <- list()
  grp_chr <- as.character(grp)
  for (i in 1:length(table(grp_chr))){
    Zgrp_index[[i]] <- which(grp_chr==names(table(ggrp_chr))[i])
  }#for i
  Zgrp <- array(0,dim=c(nrow(Z),ncol(Z),length(Zgrp_index)))
  for (i in 1:length(Zgrp_index)){
    ind <- Zgrp_index[[i]]
    Zgrp[,,i]<- Z%*%B[,ind]%*%t(B[,ind])%*%solve(B[,ind]%*%t(B[,ind]))
  }#for i
  
  Z[Z<0]=0
  Zgrp[Zgrp<0] <- 0
  
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=t(apply(Z,2,sample))
  }#if(!is.null(rseed))

  
  U=matrix(0,nrow=ncol(C), ncol=k)
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  
  iter.full.start=iter.full=20
  curfrac=0
  nposlast=Inf
  npos=-Inf
  L3.given=T
  
  for ( i in 1:max.iter){
    
    if(i>=iter.full.start){
      U=solveU(Z=Z, Chat = Chat, priorMat = C, penalty.factor = penalty.factor, glm_alpha = glm_alpha, target.frac = frac)
      curfrac=(npos<-sum(apply(U,2,max)>0))/k
      Z1=Y%*%t(B)
      for (j in 1:length(Zgrp_index)){
        ind <- Zgrp_index[[j]]
        Z1 <- Z1+Zgrp[,,j]%*%B[,ind]%*%t(B[,ind])
      }#for j
      Z2=L1*C%*%U
      ratio=median((Z2/Z1)[Z2>0&Z1>0])
      Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
    }else{
      Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    }#else
    
    Z[Z<0]=0
    
    #Zgrp
    for (j in 1:length(Zgrp_index)){
      ind <- Zgrp_index[[j]]
      Zgrp[,,j]<- Z%*%B[,ind]%*%t(B[,ind])%*%solve(B[,ind]%*%t(B[,ind]))
    }#for i
    Z[Zgrp<0] <- 0
    
    #Zgrp sparsity
    if(i>=iter.full.start & !is.null(grp)){
      for (j in 1:length(Zgrp_index)){
        Zgrp[,,j] <- BinarySearch(Zgrp[,,j], frac.grp)
      }#for j
    }#if
    
    #update B
    oldB=B
    #B=solve(t(Z-Zgrp)%*%(Z-Zgrp)+L2*diag(k))%*%t(Z-Zgrp)%*%Y
    for (j in 1:length(Zgrp_index)){
      ind <- Zgrp_index[[j]]
      B[,ind] <- solve( t(Z-Zgrp[,,j])%*%(Z-Zgrp[,,j]) + L2*diag(length(ind)))%*%t(Z-Zgrp[,,j])%*%Y[,ind]
    }#for j
    
    #doing the B sparsity
    if(i>=iter.full.start & !is.null(grp)){
      out=grpL2(tscale(B), grp)
      cutoff=quantile(out, frac.grp)
      grpMod=model.matrix(~0+grp)
      outs=scad(out, cutoff) #scad doesn't over-regularize
      outb=outs/out
      mask=outb%*%t(grpMod)
      B=B*mask
    }#if(i>=iter.full.start & !is.null(grp))
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    
    if(trace & i >=iter.full.start){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }else if (trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }#elseif
    
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }#elseif
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }#if
    
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if
    
  }#for i
  rownames(U)=colnames(priorMat)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, Zgrp=Zgrp, U=U, C=C, L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  
  if(doCrossval){
    outAUC=crossVal(out, Y, priorMat, priorMatCV)
  }else{
    message("Not using cross-validation. AUCs and p-values may be over-optimistic")
    outAUC=getAUC(out, Y, priorMat)
  }#else
  out$withPrior=which(colSums(out$U)>0)
  out$Uauc=outAUC$Uauc
  out$Up=outAUC$Upval
  out$summary=outAUC$summary
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  
  rownames(out$B)=nameB(out)
  
  out
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
}#pinv.ridge

scad=function(x, lambda,a=3.7){
  
  iip=which(abs(x)>2*lambda & abs(x)<a*lambda)
  iin=which(abs(x)<=2*lambda)
  x[iin]=pmax(0, abs(x[iin])-lambda)*sign(x[iin])
  x[iip]=((a-1)*x[iip]-sign(x[iip])*a*lambda)/(a-2)
  x
}#scad


funcByF=function (data, groups, func=mean){
  funcstats = function(x) {
    tapply(t(x), groups, func)
  }#funcstats
  t(apply(data, 1, funcstats))
}#funcByF


grpL2=function(B,grp){
  funcByF(B, grp, function(x){mean(x^2)})
}#grpl2


BinarySearch <- function (argu, grp.frac) 
{
  if (normF(argu) == 0 || sum(argu!=0)/prod(dim(argu)) <= grp.frac) 
    return(argu)
  lam1 = 0
  lam2 = max(abs(argu)) - 1e-05
  iter = 1
  while (iter < 150) {
    su = scad(argu, (lam1 + lam2)/2)
    if ( sum(argu!=0)/prod(dim(argu)) < grp.frac) {
      lam2 = (lam1 + lam2)/2
    }
    else {
      lam1 = (lam1 + lam2)/2
    }
    if ((lam2 - lam1) < 1e-01)
      #if ((lam2 - lam1) < 1e-06) 
      return(scad(argu, (lam1 + lam2)/2))
    iter = iter + 1
  }
  warning("Didn't quite converge")
  return(scad(argu, (lam1 + lam2)/2))
}



