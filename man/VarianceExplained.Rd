\name{VarianceExplained}
\alias{VarianceExplained}
\title{VarianceExplained}
\usage{
VarianceExplained(Y,Z,B, k=NULL, option = "regression")
}
\description{
Calculate explained variance for LVs decomposed by PLIER.
}
\arguments{
\item{Y}{The same data as to be used for PLIER (z-score recommended)}
\item{Z}{Loading matrix decomposed by PLIER}
\item{B}{LVs decomposed by PLIER}
\item{option}{There are four options, which are "regression", "matrix_regression", "simple" and "project". 

"project" calculates the accumulated variance and then the variance for each single LVs by subtraction. The way it calculates the accumulated variance explained is to project the input matrix onto subspace spanned by "first" several k LVs and use the trace as the variance. You can permute the LVs and corresponding loadings to see if the variance of each single LV is still consistent.

"simple", "regression" and "matrix_regression" adopt the formulation of variance explained from linear regression and they performed similarly.
}
}
\value{
A single vector with length the same as k (number of LVs) indicated by Z and B.
}
\examples{
#Y: genes-by-samples: p-by-n
#Z: genes-by-LVs: p-by-k
#B: LVs-by-samples: k-by-n
library(MASS)

#generate a toy example using SVD
P <- 2000
N <- 500
k <- 40

set.seed(1)
rnaseq <- matrix(rnorm(P*N), nrow = P)
svd.res <- svd(tscale(rnaseq), nu=k,nv=k)
Z <- svd.res$u
B <- t(svd.res$v)
PLIER.res <- list(B=B, Z=Z)

Y <- tscale(rnaseq)
Z <- PLIER.res$Z
B <- PLIER.res$B

res1 <- VarianceExplained(Y,Z,B, option="project")
plot(res1)

#Permute the LVs and corresponding loadings to see if the variance of each single LV is still consistent
#compare res2[order(order.id)] and res1
set.seed(1)
order.id <- sample(k,replace = F)
res2 <- VarianceExplained(Y,Z[,order.id],B[order.id,], option="project")
plot(res2)
plot(res2[order(order.id)])
}
