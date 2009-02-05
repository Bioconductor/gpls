\name{glpls1a.mlogit.cv.error}
\alias{glpls1a.mlogit.cv.error}
\title{Leave-one-out cross-validation error using MIRWPLS and MIRWPLSF model}
\description{
  Leave-one-out cross-validation training set error for
  fitting MIRWPLS or MIRWPLSF model for multi-group classification
}
\usage{
glpls1a.mlogit.cv.error(train.X, train.y, K.prov = NULL, eps = 0.001,lmax = 100, mlogit = T, br = T)
}

\arguments{
  \item{train.X}{ n by p design matrix (with no
    intercept term) for training set}
  \item{train.y}{  response vector with class lables 1 to C+1 for C+1 group
    classification,  baseline class should be 1}
  \item{K.prov}{ number of PLS components}
  \item{eps}{tolerance for convergence}
  \item{lmax}{ maximum number of iteration allowed }
  \item{mlogit}{if \code{TRUE} use the multinomial logit model,
    otherwise fit all C-1 logistic models (vs baseline class 1) separately}
  \item{br}{TRUE if Firth's bias reduction procedure is used}
}
\details{}

\value{
  \item{error}{LOOCV training error}
  \item{error.obs}{the misclassified error observation indices}
}

\references{
  \itemize{
    \item Ding, B.Y. and Gentleman, R. (2003) \emph{Classification using
    generalized partial least squares}. 
    \item Marx, B.D (1996) Iteratively reweighted partial least squares estimation for generalized linear regression. \emph{Technometrics} 38(4): 374-381.
  }
}

\author{Beiying Ding, Robert Gentleman}
\note{}

\seealso{ \code{\link{glpls1a.cv.error}}, \code{\link{glpls1a.train.test.error}},\code{\link{glpls1a}}, \code{\link{glpls1a.mlogit}},\code{\link{glpls1a.logit.all}}}

\examples{
 x <- matrix(rnorm(20),ncol=2)
 y <- sample(1:3,10,TRUE)

 ## no bias reduction
 glpls1a.mlogit.cv.error(x,y,br=FALSE)
 glpls1a.mlogit.cv.error(x,y,mlogit=FALSE,br=FALSE)
 ## bias reduction
 glpls1a.mlogit.cv.error(x,y,br=TRUE)
 glpls1a.mlogit.cv.error(x,y,mlogit=FALSE,br=TRUE)

}

\keyword{regression}


