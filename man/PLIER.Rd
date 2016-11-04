\name{PLIER}
\alias{PLIER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pathway-Level Information ExtractoR (PLIER): a generative
model for gene expression data
}
\description{
\deqn{\min_{U,Z,B} || Y-ZB||_F^2+\lambda_1 || Z-CU||_F^2+\lambda_2 ||B||_F^2+\lambda_3 ||U||_F^2}
\deqn{s.t. ||U||_1 \le t, rank(Z)=k, rank(B)=k}
Matrix Y would be decomposed into Z and B. Prior information C work as a regularizaiotn on the gene loading term Z. Additional sparse (L1 norm) constraint and L2 norm constraint are included in the optimization function.
}
\usage{
PLIER(data, priorMat, svdres = NULL, k, L1, L2 = NULL, L3=10, sumabsU = 5, posU = T, posZ = T, max.iter = 200, trace = F, scale = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Y, the matrix (gene expresion profile) which is going to be decomposed.
}
  \item{priorMat}{
C, the matrix records prior information.
}
  \item{svdres}{
If svdres is set to be NULL, the computational results svd(Y) would be used as the initial values for B and Z. The default value is NULL.
}
  \item{k}{
Number of latent variables.
}
  \item{L1}{
A regularization parameter.
}
  \item{L2}{
A regularization parameter which would be determined by cross-validation.
}
 \item{L3}{
A regularization parameter. 
 }
  \item{sumabsU}{
As t in the formula \eqn{||U||_1 \le t}. sumabsu controls the sparsity of U. The default value is set to be 5.
}
  \item{posU}{
If posU is set to be true, all negative entries of U will be set zero. The default value of posU is set to be true.
}
  \item{posZ}{
If posZ is set to be true, all negative entries of Z will be set zero. The default value of posZ is set to be true.
}
  \item{max.iter}{
The maximum number of iterations. The default value is set to be 200.
}
  \item{trace}{
If trace is set to be true, there will be outputs as every iteration finishes.
}
  \item{scale}{
If scale is set to true, Y will be scaled by gene.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{B}{The matrix with k latent variables}
\item{Z}{The matrix with gene loadings}
\item{C}{priorMat}
\item{U}{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
plierRes <- PLIER(data=HCCdataUse, HCCpath[cmHCC,],svdres=svdresHCC, k=52, L1=25, trace=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
