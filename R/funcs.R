ListPriors <- function(){
  output <- data(package="PLIER")[[3]]
  prior.index <- which(output[,"Title"]!="" & (!output[,"Title"] %in% c("Immune Markers (Clinical Pathways)", "Surrogate Proportion Variables","Vaccination Data" )))
  print(output[prior.index,c("Item","Title")])
}#ListPriors



VarianceExplained <- function(Y,Z,B, k=NULL, option = "regression"){
  if (is.null(k)){
    k <- min(c(ncol(Y), nrow(Y), ncol(Z), nrow(B)))
  }#if
  
  res <- c()
  
  if (option == "regression"){
    #left <- ginv(Y%*%t(Y))%*%Y
    #B and Z are required
    for (i in 1:k){
      #coef <- left%*%matrix(B[i,], ncol=1)
      coef <- matrix(Z[,i],ncol=1)
      var_reg <- sum((t(Y)%*%coef-mean(B[i,]))^2)
      res <- c(res, var_reg)
    }#for i
    
  }else if ( option == "matrix_regression"){
    #only B is required
    for (i in 1:k){
      Zreg <- Y%*%matrix(B[i,],ncol=1)%*%ginv(matrix(B[i,], nrow=1)%*%matrix(B[i,], ncol=1))
      res <- c(res, sum( ( matrix(Zreg,ncol=1)%*%matrix(B[i,],nrow=1)-mean(Y))^2))
    }#for i
  }else if (option == "simple"){
    #B and Z are required
    Z.norm <- apply(Z,2,normF)
    B.norm <- apply(B,1,normF)
    Z <- sweep(Z,2,Z.norm,"/")
    B <- sweep(B,1,B.norm, "/")
    res <- (diag(t(Z)%*%Y%*%t(B))[1:k])^2
    
  }else if (option == "project"){
    #only B is required
    res_temp <- c()
    for (i in 1:k){
      Bk <- matrix(t(B[1:i,]), ncol=i)
      Xk <- Y%*%Bk%*%ginv(t(Bk)%*%Bk)%*%t(Bk)
      res_temp <- c(res_temp, sum(diag(t(Xk)%*%Xk)))  
    }#for i
    res <- res_temp
    res[2:k] <- res[2:k]-res_temp[1:(k-1)]
    
  }#else if
  return(res)
}#VarianceExplained
