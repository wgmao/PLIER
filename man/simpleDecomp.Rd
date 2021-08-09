\name{simpleDecomp}
\alias{simpleDecomp}
\title{PLIER without prior}
\usage{
 simpleDecomp=function(Y, k,svdres=NULL, L1=NULL, L2=NULL,
         Zpos=T,max.iter=200, tol=5e-6, trace=F,
         rseed=NULL, B=NULL, scale=1, pos.adj=3, cutoff=0)
}
\arguments{
  \item{Y}{the data to be processed with genes in rows and samples in columns. Should be z-scored or set scale=T}
  
  \item{k}{The number of latent variables to return}
  
  \item{svdres}{Pre-computed result of the svd decomposition for data}
  
  \item{L1}{L1 constant, leave as NULL to automatically select a value}
  
  \item{L2}{L2 constant, leave as NULL to automatically select a value}
  
  \item{Zpos}{If Zpos is set to be T, values in the loading matrix which are smaller than cutoff will be set to zeroes. }
  
  \item{max.iter}{Maximum number of iterations to perform}
  
  \item{tol}{Convergence threshold}
  
  \item{trace}{Display progress information}
  
  \item{rseed}{Set this option to use a random initialization, instead of SVD}
  
  \item{B}{Initial values for B (latent variables)}
  
  \item{scale}{Coefficient to adjust L1 and L2}
  
  \item{pos.adj}{Coefficient to adjust L1}
  
  \item{cutoff}{see Zpos}
  
\description{
  PLIER decomposition without prior
}