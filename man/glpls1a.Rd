\name{glpls1a}
\alias{glpls1a}
\title{Fit IRWPLS and IRWPLSF model}
\description{
  Fit Iteratively ReWeighted Least Squares (IRWPLS) with an option of
  Firth's bias reduction procedure (IRWPLSF) for two-group classification
}
\usage{
glpls1a(X, y, K.prov = NULL, eps = 0.001, lmax = 100, b.ini = NULL, 
      denom.eps = 1e-20, family = "binomial", link = NULL, br = TRUE)
}

\arguments{
  \item{X}{ n by p design matrix (with no intercept term)}
  \item{y}{ response vector 0 or 1}
  \item{K.prov}{ number of PLS components, default is the rank of X}
  \item{eps}{tolerance for convergence}
  \item{lmax}{ maximum number of iteration allowed }
  \item{b.ini}{ initial value of regression coefficients}
  \item{denom.eps}{ small quanitity to guarantee nonzero denominator in
    deciding convergence}
  \item{family}{ glm family, \code{binomial} is the only relevant one here }
  \item{link}{ link function, \code{logit} is the only one practically implemented now}
  \item{br}{TRUE if Firth's bias reduction procedure is used}
}
\details{
  
}
\value{
  \item{coefficients }{regression coefficients}
  \item{convergence }{whether convergence is achieved}
  \item{niter}{total number of iterations}
  \item{bias.reduction}{whether Firth's procedure is used}
}

\references{
 \item Ding, B.Y. and Gentleman, R. (2003) Classification using generalized partial least squares.
 \item Marx, B.D (1996) Iteratively reweighted partial least squares estimation for generalized linear regression. Technometrics 38(4): 374-381.
}
\author{Beiying Ding, Robert Gentleman}
\note{}

\seealso{ \code{\link{glpls1a.mlogit}}, \code{\link{glpls1a.logit.all}},
  \code{\link{glpls1a.train.test.error}},
  \code{\link{glpls1a.cv.error}}, \code{\link{glpls1a.mlogit.cv.error}}}

\examples{
 x <- matrix(rnorm(20),ncol=2)
 y <- sample(0:1,10,TRUE)
 ## no bias reduction
 glpls1a(x,y,br=FALSE)
  
 ## no bias reduction and 1 PLS component
 glpls1a(x,y,K.prov=1,br=FALSE)

 ## bias reduction
 glpls1a(x,y,br=TRUE)
}

\keyword{regression}


