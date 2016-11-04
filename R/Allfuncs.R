library(qvalue)
library(RColorBrewer)
library(pheatmap)
library(survival)

rowNorm=function(x){
  s=apply(x,1,sd)
  m=apply(x,1,mean);
  x=sweep(x,1,m)
  x=sweep(x,1,s,"/")
  x
}

PLIER=function(data, priorMat,svdres=NULL, k, L1, L2=NULL, L3=10, sumabsU=5, posU=T, posZ=T, max.iter=200, trace=F,   scale=F){
  if(is.null(L2)){
    L2=nrow(data)/k
    print(paste0("L2 is set to ",L2))
  }
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
  
  Chat=pinv.ridge(crossprod(C), L3)%*%(t(C))
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  if(is.null(svdres)){
    message("Computing SVD")
    svdres=svd(Y)
    message("Done")
  }
  B=t(svdres$v[, 1:k]%*%diag(svdres$d[1:k]))
  Z=svdres$u[, 1:k]
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  
  for ( i in 1:max.iter){
    
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace){
      
      message(paste0("errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", errorZ= ",errz<-round2(mean((Z-C%*%U)^2)), ", ratio= ", round(erry/errz), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    U=Chat%*%Z
    if(posU){
      U[U<0]=0
    }
    
    U=BinarySearch(U,sumabsU)
    
    Z=(Y%*%t(B)+L1*C%*%U)%*%solve(tcrossprod(B)+L1*diag(k))
    
    if(posZ){
      Z[Z<0]=0
    } 
    
    oldB=B
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
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
  return(list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C))
}




testSurvival=function(features,ctime, cstatus, covariates=NULL){
  ii=which(cstatus!="" & !is.na(ctime))
  cdata=list(time=ctime[ii], status=as.numeric(as.factor(as.character(cstatus[ii]))))
  out=matrix(nrow=nrow(features), ncol=2)
  if (!is.null(covariates)){
    for( i in 1:nrow(features)){
      cdata$x=t(covariates[, ii, drop=F])
      cfit0<-coxph(Surv(time, status) ~ x, cdata, model=T);
      cdata$x=cbind(cdata$x,features[i,ii])
      cfit<-coxph(Surv(time, status) ~ x, cdata, model=T);
      res=anova(cfit,cfit0)
      out[i,1]=res$`P(>|Chi|)`[2]
      out[i,2]=cfit$coefficients[nrow(covariates)+1]
    }
  }
  else{
    for( i in 1:nrow(features)){
      cdata$x=features[i,ii]
      cfit<-coxph(Surv(time, status) ~ x, cdata, model=T);out[i,1]=(summary(cfit)$logtest[3])
      out[i,2]=cfit$coefficients
      
    }
  }
  QV= function(x){(qvalue(x))$qvalue}
  list(p.val=out[,1], coef=out[,2], q.value=QV(out[,1]))
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




combineGSmats=function(...){
  GSmats=list(...)
  genes=character()
  
  for( i in 1:length(GSmats)){
    genes=c(genes, rownames(GSmats[[i]]))
    GSmats[[i]]=as.data.frame(GSmats[[i]])
  }
  
  genes=unique(genes)

  matAll=matrix(nrow=length(genes), ncol=0)
  rownames(matAll)=genes
  for( i in 1:length(GSmats)){
    GSmats[[i]]=as.matrix(GSmats[[i]][genes,])
    matAll=cbind(matAll, GSmats[[i]])
  }
  matAll[is.na(matAll)]=0
  
  matAll
}


commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
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


l2n = function (vec) 
{
  a <- sqrt(sum(vec^2))
  if (a == 0) 
    a <- 0.05
  return(a)
}


soft = function (x, d) 
{
  return(sign(x) * pmax(0, abs(x) - d))
}


BinarySearch=function (argu, sumabs) 
{
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
    return(argu)
  lam1 = 0
  lam2 = max(abs(argu)) - 1e-05
  iter = 1
  while (iter < 150) {
    su = soft(argu, (lam1 + lam2)/2)
    if (sum(abs(su/l2n(su))) < sumabs) {
      lam2 = (lam1 + lam2)/2
    }
    else {
      lam1 = (lam1 + lam2)/2
    }
    if ((lam2 - lam1) < 1e-06) 
      return(soft(argu, (lam1 + lam2)/2))
    iter = iter + 1
  }
  warning("Didn't quite converge")
  return(soft(argu, (lam1 + lam2)/2))
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