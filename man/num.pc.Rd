\name{num.pc}
\alias{num.pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count the number of latent variables in the gene expression profile
}
\description{
Count the number of latent variables in the gene expression profile which could help to determin the parameter k in the main function PLIER().
}
\usage{
num.pc(data, B = 20, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The gene expression profile matrix, or matrix Y in PLIER().
}
  \item{B}{
The number of times dat would be resampled.
}
  \item{seed}{
seed is used to set the random number generator. The default value is set to be NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The number of possible latent variables embedded in the gene expression profile.
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
res <- num.pc(HCCdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
