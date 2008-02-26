\name{predict.gpls}
\alias{predict.gpls}
\alias{predict,gpls-method}
\title{A prediction method for gpls.}
\description{
  A simple prediction method for \code{gpls} objects.
}
\usage{
\method{predict}{gpls}(object, newdata, ...)
}
\arguments{
  \item{object}{A \code{gpls} object, typically obtained from a call to
  \code{\link{gpls}}  }
  \item{newdata}{New data, for which predictions are desired. }
  \item{\dots}{Other arguments to be passed on}
}
\details{
  The prediction method is straight forward. The estimated coefficients
  from \code{object} are used, together with the new data to produce
  predicted values. These are then split, according to whether the
  predicted values is larger or smaller than 0.5 and predictions
  returned.

  The code is similar to that in \code{\link{glpls1a.train.test.error}}
  except that in that function both the test and train matrices are
  centered and scaled (the covariates) by the same values (those from
  the test data set).
}
\value{
  A \code{list} of length two:
  \item{class}{The predicted classes; one for each row of \code{newdata}.}
  \item{predicted}{The estimated predictors.}
}
\author{B. Ding and R. Gentleman}
\seealso{\code{\link{gpls}}}
\examples{
  example(gpls)
  p1 = predict(m1)
}
\keyword{classif}
