\name{gpls}
\alias{gpls}
\alias{print.gpls}
\alias{gpls.formula}
\alias{gpls.default}

\title{A function to fit Generalized partial least squares models. }
\description{
 Partial least squares is a commonly used dimension reduction
 technique. The paradigm can be extended to include generalized linear
 models in several different ways. The code in this function uses the
 extension proposed by Ding and Gentleman, 2004.
}
\usage{
gpls(x, \dots)

\method{gpls}{default}(x, y, K.prov=NULL, eps=1e-3, lmax=100, b.ini=NULL,
    denom.eps=1e-20, family="binomial", link=NULL, br=TRUE, \dots)

\method{gpls}{formula}(formula, data, contrasts=NULL, K.prov=NULL,
eps=1e-3, lmax=100, b.ini=NULL, denom.eps=1e-20, family="binomial",
link=NULL, br=TRUE, \dots)

}

\arguments{
  \item{x}{The matrix of covariates.}
  \item{formula}{ A formula of the form 'y ~ x1 + x2 + ...', where
    \code{y} is the response and the other terms are covariates.}
  \item{y}{The vector of responses}
  \item{data}{A data.frame to resolve the forumla, if used}
  \item{K.prov}{ number of PLS components, default is the rank of X}
  \item{eps}{tolerance for convergence}
  \item{lmax}{ maximum number of iteration allowed }
  \item{b.ini}{ initial value of regression coefficients}
  \item{denom.eps}{ small quanitity to guarantee nonzero denominator in
    deciding convergence}
  \item{family}{ glm family, \code{binomial} is the only relevant one here }
  \item{link}{ link function, \code{logit} is the only one practically
    implemented now} 
  \item{br}{TRUE if Firth's bias reduction procedure is used}
  \item{\dots}{Additional arguements. }
  \item{contrasts}{ an optional list. See the \code{contrasts.arg} of
          \code{model.matrix.default}.}
}
\details{
  This is a different interface to the functionality provided by
  \code{\link{glpls1a}}. The interface is intended to be simpler to use
    and more consistent with other matchine learning code in R.

    The technology is intended to deal with two class problems where
    there are more predictors than cases. If a response variable
    (\code{y}) is used that has more than two levels the behavior may
    be unusual.
}
\value{
  An object of class \code{gpls} with the following components:
  \item{coefficients}{The estimated coefficients.}
  \item{convergence}{A boolean indicating whether convergence was
    achieved.}
  \item{niter}{The total number of iterations.}
  \item{bias.reduction}{A boolean indicating whether Firth's procedure
    was used.} 
  \item{family}{The \code{family} argument that was passed in.}
  \item{link}{The \code{link} argument that was passed in.}
  \item {terms}{The constructed terms object.}
  \item{call}{The call}
  \item{levs}{The factor levels for prediction.}
}
\references{  \item{Ding, B.Y. and Gentleman, R. (2003)}{Classification
using generalized partial least squares.}
 \item{Marx, B.D (1996)}{Iteratively reweighted partial least squares
 estimation for generalized linear regression. Technometrics 38(4):
 374-381.}
}
\author{B. Ding and R. Gentleman }

\seealso{ \code{\link{glpls1a}} }
\examples{
library(MASS)
m1 = gpls(type~., data=Pima.tr, K=3)

}
\keyword{ classif }

