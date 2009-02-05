\name{glpls1a.train.test.error}
\alias{glpls1a.train.test.error}
\title{out-of-sample test set error using IRWPLS and IRWPLSF model}
\description{
  Out-of-sample test set error for fitting IRWPLS or IRWPLSF model on the
  training set for two-group classification
}
\usage{
glpls1a.train.test.error(train.X,train.y,test.X,test.y,K.prov=NULL,eps=1e-3,lmax=100,family="binomial",link="logit",br=T)
}

\arguments{
  \item{train.X}{  n by p design matrix (with no
    intercept term) for training set}
  \item{train.y}{ response vector (0 or 1) for training set}
  \item{test.X}{ transpose of the design matrix (with no
    intercept term) for test set}
  \item{test.y}{ response vector (0 or 1) for test set}
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
  \item{error}{out-of-sample test error}
  \item{error.obs}{the misclassified error observation indices}
  \item{predict.test}{the predicted probabilities for test set}
}

\references{
  \itemize{
    \item Ding, B.Y. and Gentleman, R. (2003) \emph{Classification using
    generalized partial least squares}. 
    \item Marx, B.D (1996) Iteratively reweighted partial least squares
    estimation for generalized linear regression. \emph{Technometrics} 38(4):
    374-381. 
  }
}
\author{Beiying Ding, Robert Gentleman}
\note{}

\seealso{ \code{\link{glpls1a.cv.error}},
  \code{\link{glpls1a.mlogit.cv.error}}, \code{\link{glpls1a}}, \code{\link{glpls1a.mlogit}}, \code{\link{glpls1a.logit.all}}}

\examples{
 x <- matrix(rnorm(20),ncol=2)
 y <- sample(0:1,10,TRUE)
 x1 <- matrix(rnorm(10),ncol=2)
 y1 <- sample(0:1,5,TRUE)

 ## no bias reduction
 glpls1a.train.test.error(x,y,x1,y1,br=FALSE)
 ## bias reduction
 glpls1a.train.test.error(x,y,x1,y1,br=TRUE)
}

\keyword{regression}


