glpls1a.logit.all <- function(X,y,K.prov=NULL,eps=1e-3,lmax=100,b.ini=NULL,denom.eps=1e-20,family="binomial",link="logit",br=T)
{
  family <-"binomial"
  link <- "logit"

  x <- as.matrix(X)
  yv <- as.numeric(as.vector(y))
  dimnames(x) <- names(yv) <- NULL
  C <- max(yv)-1

  n <- nrow(x)
  r <- ncol(x)+1
  
  y <- matrix(0,n,C+1)
  y[cbind(seq(n),yv)] <- 1
  y0 <- y[,1]
  y <- y[,-1,drop=F]

  beta <- matrix(0,r,C)

  for ( i in 1:C)
    {
      index <- (1:n)[y0 != y[,i]]
      X <- x[index,]
      Y <- as.vector(y[index,i])
      beta[,i] <- glpls1a(X,Y,K.prov=K.prov,eps=eps,lmax=lmax,b.ini=b.ini,denom.eps=denom.eps,family=family,link=link,br=br)$coef
    }
  return(list(coefficients=beta))
}
