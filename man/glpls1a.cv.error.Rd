\name{glpls1a.cv.error}
\alias{glpls1a.cv.error}
\title{Leave-one-out cross-validation error using IRWPLS and IRWPLSF model}
\description{
  Leave-one-out cross-validation training set classification error for
  fitting IRWPLS or IRWPLSF model for two group classification
}
\usage{
glpls1a.cv.error(train.X,train.y, K.prov=NULL,eps=1e-3,lmax=100,family="binomial",link="logit",br=T)
}

\arguments{
  \item{train.X}{ n by p design matrix (with no
    intercept term) for training set}
  \item{train.y}{ response vector (0 or 1) for training set}
  \item{K.prov}{ number of PLS components, default is the rank of
    train.X}
  \item{eps}{tolerance for convergence}
  \item{lmax}{ maximum number of iteration allowed }
  \item{family}{ glm family, \code{binomial} is the only relevant one here }
  \item{link}{ link function, \code{logit} is the only one practically implemented now}
  \item{br}{TRUE if Firth's bias reduction procedure is used}
}
\details{
  
}
\value{
  \item{error}{LOOCV training error}
  \item{error.obs}{the misclassified error observation indices}
}

\references{
 \item Ding, B.Y. and Gentleman, R. (2003) Classification using generalized partial least squares.
 \item Marx, B.D (1996) Iteratively reweighted partial least squares estimation for generalized linear regression. Technometrics 38(4): 374-381.
}
\author{Beiying Ding, Robert Gentleman}
\note{}

\seealso{ \code{\link{glpls1a.train.test.error}},
  \code{\link{glpls1a.mlogit.cv.error}}, \code{\link{glpls1a}}, \code{\link{glpls1a.mlogit}},\code{\link{glpls1a.logit.all}}}

\examples{
 x <- matrix(rnorm(20),ncol=2)
 y <- sample(0:1,10,TRUE)

 ## no bias reduction
 glpls1a.cv.error(x,y,br=FALSE)
 ## bias reduction and 1 PLS component
 glpls1a.cv.error(x,y,K.prov=1, br=TRUE)
}

\keyword{regression}


