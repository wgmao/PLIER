require(RColorBrewer)
require(gplots)
require(pheatmap)
require(glmnet)
require(rsvd)
require(qvalue)

#' @keywords  internal
QV<-function(pval){
  
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


rsvd=function(...){
  message("Using rsvd: set seed for consistent results")
rsvd::rsvd(...)
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
  iimulit=which(iin>1)
 
  message(paste0("Genes not matched uniquely: ", paste(names(iimulit), collapse=", ")))

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

#' names latent variables according to their pathway useage (if any)
#' 
#' @param plierRes A PLIER result
#' @param top The number of pathway to use. Only the top pathway (one with the largest coefficient) is used by default
#' @param fdr.cutoff The cross-validation significance cutoff for a pathway to be considered for naming
#' @export
nameB=function(plierRes, top=1, fdr.cutoff=0.05){

  names=vector("character",ncol(plierRes$U))
  Uuse=plierRes$U
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
    
    for(j in iipath){
      iiheldout=which((priorMat[,j]==0) |(priorMat[,j]==1&priorMatcv[,j]==0))
      aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
      out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=QV(out[,4])
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
#' @param computeAUC Whether or not to do pseudo cross-validation
#' @param penalty.factor A vector equal to the number of columns in priorMat. Sets relative penalties for different pathway/geneset subsets. Lower penalties will make a pathway more likely to be used. Only the relative values matter. Internally rescaled.
#' @param glm_alpha Set the alpha for elastic-net
#' @param minGenes The minimum number of genes a pathway must have to be considered
#' @param tol Convergence threshold
#' @param seed Set the seed for reproducible results
#' @export

PLIER=function(data, priorMat,svdres=NULL, k=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=F, Chat=NULL, maxPath=10, computeAUC=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-5, seed=NULL){
  
  #Ur is the ranked matrix of pathway relevance
  solveU=function(Z, Ur,priorMat,  L3, penalty.factor, glm_alpha){
    
  
    ii=which(apply(Ur,1,min)<=maxPath)
    U=copyMat(Ur)
    U[]=0
    for (j in 1:ncol(U)){
      selection=which(Ur[,j]<=maxPath)
      tmp=glmnet(y=Z[,j], x=priorMat[,selection], alpha=glm_alpha, lambda=L3, lower.limits = 0, penalty.factor = penalty.factor[selection])
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
  
  if(!all(rownames(priorMat)==rownames(data))){
    cm=commonRows(data, priorMat)
    warning(paste("Selecting common genes:", length(cm)))
    priorMat=priorMat[cm,]
    data=data[cm,]

  }
  numGenes=colSums(priorMat)
  iibad=which(numGenes<minGenes)
  warning(paste("Removing", length(iibad), "pathways with too few genes"))
  priorMat[, iibad]=0
  heldOutGenes=list()
if(computeAUC){
  priorMatCV=priorMat
  if(!is.null(seed))
    set.seed(seed)
  for(j in 1:ncol(priorMatCV)){
    iipos=which(priorMatCV[,j]==1)
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
    Chat=pinv.ridge(Cp, 5)%*%(t(C))
  }
  YsqSum=sum(Y^2)
  #compute svd and use that as the starting point

  
  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
    svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }
    else{
      svdres=svd(Y)
    }
    message("Done")
  }
  if(is.null(k)){
    k=num.pc(svdres)
    message("k is set to ", k)
  }
  if(nrow(svdres$u)!=nrow(data)){
    message("SVD U has the wrong number of rows")
    
    if(!is.null(rownames(svdres$u))){
      message("Selecting via rownames")
      Z=svdres$u[rownames(Y),1:k]
    }
    else{
      message("No rownames for svdres$u: Recomputing SVD")
      if(ns>500){
        if(!is.null(seed))
          set.seed(seed)
        svdres=rsvd(Y, k, q=3)
      }
      else{
        svdres=svd(Y)
      }
      
      message("Done")
    }
  }

    Z=svdres$u[, 1:k]
  
 
  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }
  
  
  
  B=t(svdres$v[, 1:k]%*%diag(svdres$d[1:k]))
  U=matrix(0,nrow=ncol(C), ncol=k)


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
      ii=which(apply(Us,1,min)<=maxPath)
      
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
          
          U=solveU(Z, Us, priorMat,  L3=(L3use<-(L3_1+L3_2)/2), penalty.factor, glm_alpha)
          
          nposlast=npos
          curfrac=(npos<-sum(apply(U,2,max)>0))/k
         # message(paste0(npos, " positive columns at L3=", round(L3use,6)))
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
      ii=which(apply(Us,1,min)<=maxPath)
      U=solveU(Z, Us, priorMat, L3, penalty.factor, glm_alpha)
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

  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, numActPath=length(ii), L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)

  if(computeAUC){
    crossval=crossVal(out, Y, priorMat, priorMatCV)
    
    out$withPrior=which(colSums(out$U)>0)
    out$Uauc=crossval$Uauc
    out$Up=crossval$Upval
    out$summary=crossval$summary
    tt=apply(out$Uauc,2,max)
    message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  }
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
#' @params ... Additional arguments to be passed to pheatmap, such as a column annotation data.frame
#' @export
plotTopZ=function(plierRes, data, priorMat, top=10, index=NULL, regress=F, allLVs=F,...){
  
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
      toPlot[gi,]=toPlot[gi, ]-plierRes$Z[rownames(toPlot)[gi],-i ]%*%plierRes$B[-i,colnames(toPlot)]
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
  if(class(data)!="list"){
  n <- ncol(data)
  m <- nrow(data)
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
  uu <- rsvd(data,k, q=3)
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
      uu <- svd(data)
    }
    else{
      uu <- rsvd(data,k, q=3)
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
  }
  else if (method=="elbow"){
    x<-smooth(abs(diff(diff(uu$d))))
    nsv=which(x<median(x))[1]
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
  
  mydist=function(x){dist(abs(sign(x)))}
  
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
 
  pheatmap(matrix,color =mycol , clustering_callback = function(h,d){hclust(mydist(d), method = "single")}, ...)

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
#' @params ... Additional arguments to be passed to pheatmap, such as a column annotation data.frame
#' @export
plotTopZallPath=function(plierRes, data, priorMat, top=10, index=NULL, regress=F,fdr.cutoff=0.05,...){
  
  
  pval.cutoff=max(plierRes$summary[plierRes$summary[,5]<fdr.cutoff,4])
  ii=which(colSums(plierRes$U)>0)
  if(! is.null(index)){
    ii=intersect(ii,index)
  }
  tmp=apply(-plierRes$Z[, ii, drop=F],2,rank)
  nn=character()
  nncol=character()
  nnpath=character()
  nnindex=double()
  Ustrict=plierRes$U
  Ustrict[plierRes$Up>pval.cutoff]=0
  
  pathsUsed=which(rowSums(Ustrict[, index, drop=F])>0)
  pathMat=matrix(nrow=0, ncol=length(pathsUsed))
 
  colnames(pathMat)=strtrim(names(pathsUsed),30)

  for (i in 1:length(ii)){
    nn=c(nn,nntmp<-names(which(tmp[,i]<=top)))
    nncol=c(nncol, rep(rownames(plierRes$U)[which(thispath<-plierRes$U[,ii[i]]==max(plierRes$U[,ii[i]]))], length(nntmp)))
    nnindex=c(nnindex,rep(ii[i], length(nntmp)))
    
    pathMat=rbind(pathMat, priorMat[nntmp, pathsUsed])
  }
  pathMat=pathMat[, colSums(pathMat)>0]
  pathMat=as.data.frame(pathMat)
  pathMat=apply(pathMat,2,as.factor)
  names(nncol)=nn
  nncol=strtrim(nncol, 30)
  
  nnrep=names(which(table(nn)>1))
  

  ll=list(inPathway="black", notInPathway="beige")
 
  ll2=list()
  for(i in 1:length(pathsUsed)){
    ll2[[i]]=c("black", "beige")
    names(ll2[[i]])=c("1","0")
  }
  
  names(ll2)=colnames(pathMat)
  

  anncol=ll2

  mydist=function(x){as.dist(1-t(cor(t(x))))}
  rr=max(range(tscale(data[nn,])))
  bb=seq(-rr,rr, length.out = 100)

  toPlot=data[nn,]
  if(regress){
    for ( i in ii){
      gi=which(nnindex==i)
  
      toPlot[gi,]=toPlot[gi, ]-plierRes$Z[rownames(toPlot)[gi],-i ]%*%plierRes$B[-i,colnames(toPlot)]
    }
  }


#  show(pathMat)
#  pheatmap(tscale(toPlot), breaks=bb,color=colorpanel(101, "green", "white", "red"),annotation_row=as.data.frame(pathMat))
  pheatmap(tscale(toPlot), breaks=bb,color=colorpanel(101, "green", "white", "red"),annotation_row=as.data.frame(pathMat[, ]), annotation_legend = F, 
           show_colnames = F, annotation_colors = anncol,
           clustering_callback = function(h,d){hclust(mydist(d), method = "average")},...)
}
