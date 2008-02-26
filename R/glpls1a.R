##Copyright 2004 B. Ding and R. Gentleman, all rights reserved

glpls1a <- function(X, y, K.prov=NULL, eps=1e-3, lmax=100, b.ini=NULL,
                    denom.eps=1e-20, family="binomial", link=NULL, br=TRUE)
{

  if (is.null(link)) {
    if (family=="normal") link <- "identity"
    else if (family=="binomial") link <- "logit"
    else if (family=="poisson") link <- "log"
    else stop("unknown family", family)
  }

  X <- as.matrix(X)
  dx <- dim(X)

##FIX up Y
  if (is.factor(y)) {
      levs = levels(y)
      if( length(levs) > 2 )
          levs = c(levs[1], "Other")
      y <- y != levels(y)[1]
  }
  else levs = unique(y)
  if(length(levs) > 2)
      warning("y has more than two levels")

  if (any(y < 0 | y > 1))
        stop("y values must be 0 <= y <= 1")

### number of PLS components

  if (is.null(K.prov)) K <- min(dx[1]-1,dx[2])
  else if ( K.prov > min(dx[1]-1,dx[2]) )
    {
      cat("number of dimension K.provd exceeds rank of X. \n","Provide a number that is less than or equal to ",min(dx[1]-1,dx[2]),"!\n");
      K <- min(dx[1]-1,dx[2])
    }
  else
    {
      K <- K.prov
    }

  ## cat("Number of components is:", K,"!\n")

### initialzing matrices for PLS regression

  ## weight matrix for PLS regression
  W<- matrix(0,dx[2],K)

  ## score matrix
  Ta <- matrix(0,dx[1],K)

  ## loading matrix for X
  P <- matrix(0,dx[2],K)

  ## loading matrix for y
  Q <- numeric(K)

### number of iterations
  l <- 0

### intial values

  ## convergence
  converged <- F

  ## initial predictor matrix
  E <- X

  ## pseudo response
  ystar <- psi(y, family)
  ystar0 <- ystar

  ## Weight matrix for GLM, e.g. p(1-p) for logistic regression
  V <- diag(c(hp(ystar,family,link)^2/bpp(ystar,family,link)))

  ## weighted standardization of predictor matrix
  E <- t(E)-apply(E,2,weighted.mean,diag(V))
  E <- t(E)
#  E <- t(E/sqrt(apply(E,1,var)))

  ## linear predictor, ie. eta = Xb
  eta <- rep(0,length(y))
  eta.old <- eta+10

  Q <- rep(0,K)
  Q.old <- rep(100,K)
  K.old <- K
  min <- 1000

  ## regression coefficients
  if(!is.null(b.ini))
    {
      beta.old <- b.ini[-1]
    }
  else{
    beta <- rep(1,ncol(X))
  }

  beta.old <- beta/1000

  while ((max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)) > eps) &
         (l < lmax))     {
     ##    cat("Iter:",l,"\n")

      if ((max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)) < min))
        {
          if(l==1)
            {
              W.min <- W
              P.min <- P
             Q.min <- Q
             }
          else
            {
              ## record minimum values
              min <- max(abs(beta-beta.old)/(abs(beta.old)+denom.eps))
              W.min <- W
              P.min <- P
              Q.min <- Q
              eta.min <- eta
              l.min <- l
            }
        }

      K.old <- K

      l <- l + 1

   #   K <- min(dx[1]-1,dx[2],Rank(E),K.prov)
      W <- matrix(0,dx[2],K)
      Ta <- matrix(0,dx[1],K)
      P <- matrix(0,dx[2],K)
      Q <- numeric(K)

### PLS regression

      for ( i in 1:K) {
        w <- t(E)%*%V%*%ystar
        w <- w/sqrt(crossprod(w)[1])
        W[,i] <- w
        ta <- E%*%w
        ## ta <- ta-weighted.mean(ta,diag(V))
        ## ta <- ta/sqrt(var(ta))[1]
        Ta[,i] <- ta
        taa <- (t(ta)%*%V%*%ta)[1]
        p <- t(E)%*%V%*%ta/taa
        P[,i] <- p
        q <- (t(ystar)%*%V%*%ta/taa)[1]
        Q[i] <- q
        ystar <- ystar - q*ta
        E <- E - ta%*%t(p)
      }

      ## update beta
      beta.old <- beta
      beta <- W%*%solve(t(P)%*%W)%*%Q

      ## Hat matrix
      H <- hat(sweep(cbind(rep(1,dx[1]),X),1,sqrt(diag(V)),"*"),intercept=FALSE)
      ## update eta (linear predictor)
      eta <- weighted.mean(ystar0,diag(V)) + Ta%*%Q

      ## update and rescaling weight matrix
      V <- diag(c(hp(eta,family,link)^2/bpp(eta,family,link)))
      V <- V*(H*br+1)

      ## diagnosis for divergence
      if( sum(diag(V)=="NaN") > 0)
        {
          #print("diagonal elements of V overflow!")
           break
        }
      if (sum(round(probcal(y,eta),4)>=0.9999)==length(y))
        {
           #print("complete separation !")
           break
        }

      ## update pseudo response
     ystar <- eta +
     diag(c(1/hp(eta,family,link)))%*%(y+H*br/2 -
     (H*br+1)*h(eta,family,link))/(H*br +1)
      ystar0 <- ystar

      ## update predictor matrix
      E <- t(X)-apply(X,2,weighted.mean,diag(V))
      E <- t(E)
      ## E <- t(E/sqrt(apply(E,1,var)))

#         print(max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)))
    }

#  print(l)

  if (max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)) > eps)
    {
 #     cat("Convergence Not achieved and estimate from iteration ", l.min, " is used!")
        W <- W.min
        P <- P.min
        Q <- Q.min
        eta <- eta.min
    }
  else
    {
      converged <- T
    }

  ## final estimates
  beta <- W%*%solve(t(P)%*%W)%*%Q
  beta0 <- eta[1]-X[1,]%*%beta

  ##put some names onto the coefs
  coef = c(beta0,beta)
  dnx = dimnames(X)[[2]]
  if(length(dnx) != length(beta))
      dnx = paste("X", 1:length(beta), sep=":")
  names(coef) = c("Intercept", dnx)

  ans = list(coefficients=coef,
                convergence = converged,
                niter = l,
                family = family,
                link = link,
                levs = levs,
                bias.reduction = br)
  class(ans) = "gpls"
  ans
}

##let's see if we can put this into a standard modeling framework
gpls = function(x, ...)
    UseMethod("gpls")

gpls.default = function(x, y, K.prov=NULL, eps=1e-3, lmax=100, b.ini=NULL,
denom.eps=1e-20, family="binomial", link=NULL, br=TRUE, ...)
       glpls1a(x, y, K.prov, eps, lmax, b.ini, denom.eps, family,
       link, br)

gpls.formula = function(formula, data, contrasts=NULL, K.prov=NULL,
eps=1e-3, lmax=100, b.ini=NULL, denom.eps=1e-20, family="binomial",
link=NULL, br=TRUE, ...) {
    mf = match.call()
    m = match(c("formula", "data"), names(mf), 0)
    mf = mf[c(1,m)]
    mf[[1]] = as.name("model.frame")
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    y = model.response(mf, "numeric")
    x = if( !is.empty.model(mt))
           model.matrix(mt, mf, contrasts)
        else matrix(, NROW(y), 0)
    xint = match("(Intercept)", colnames(x), nomatch=0 )
    if(xint > 0 )
        x <- x[, -xint, drop=FALSE]
    ans = glpls1a(x, y, K.prov, eps, lmax, b.ini,
            denom.eps, family, link, br)
    ans$terms = mt
    ans$call = match.call()
    ans
}

print.gpls = function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
            quote = FALSE)
    }
    invisible(x)
}

##Based on predict.lda
predict.gpls = function (object, newdata, ...)
{
    if (!inherits(object, "gpls"))
        stop("object not of class gpls")
    if (!is.null(Terms <- object$terms)) {
        Terms <- delete.response(Terms)
        if (missing(newdata))
            newdata <- model.frame(object)
        else {
            newdata <- model.frame(Terms, newdata, na.action = na.pass,
                xlev = object$xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses"))) 
                stats::.checkMFClasses(cl, newdata)
        }
        x <- model.matrix(Terms, newdata, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0)
        if (xint > 0)
            x <- x[, -xint, drop = FALSE]
    }
    else {
        if (missing(newdata)) {
            if (!is.null(sub <- object$call$subset))
                newdata <- eval.parent(parse(text =
                paste(deparse(object$call$x, backtick = TRUE), "[",
                deparse(sub, backtick = TRUE), ",]")))
            else newdata <- eval.parent(object$call$x)
            if (!is.null(nas <- object$call$na.action))
                newdata <- eval(call(nas, newdata))
        }
        if (is.null(dim(newdata)))
            dim(newdata) <- c(1, length(newdata))
        x <- as.matrix(newdata)
    }
    if (ncol(x) != length(object$coefficients) - 1 )
        stop("wrong number of variables")
    if (length(colnames(x)) > 0 && any(colnames(x) !=
                names(object$coef[-1])) )
        warning("Variable names in newdata do not match those in object")
##FIXME: just copied this from Beiying's stuff - clearly needs a bit
## more thought - why do we strip off the coef when fitting, and then
## add it back in?
    eta <- cbind(rep(1,nrow(x)),x)%*%object$coefficient
    preds = h(eta, object$family, object$link)
    ##we seem to need to do a lot of futzing to get the right levels
    rf = preds > 0.5
    if(all(rf) ) {
        rf = factor(rf, labels=object$levs[2])
        levels(rf) = object$levs
    }
    else if( !any(rf) ) {
        rf = factor(rf, labels=object$levs[1])
        levels(rf) = object$levs
    }
    else
        rf = factor(preds[,1]>0.5, labels=object$levs)
    ans = list(class = rf, predicted = preds)
    ans
}
