glpls1a.mlogit <- function(x,y,K.prov=NULL,eps=1e-3,lmax=100,b.ini=NULL,denom.eps=1e-20,family="binomial",link="logit",br=T)
{

  if(family != "binomial" | link != "logit")
    {
      print("wrong family (link) !\n")
      break
    }

  if(any(x[,1] !=1))
    x <- cbind(rep(1,nrow(x)),x)
  
  x <- as.matrix(x)
  yv <- as.numeric(as.vector(y))
  dimnames(x) <- names(yv) <- NULL
  r <- ncol(x)
  C <- max(yv)-1

  n <- nrow(x)
  y <- matrix(0,n,C+1)
  y[cbind(seq(n),yv)] <- 1
  y0 <- y[,1]
  y <- y[,-1]
  jn <- rep(1,n)
  jC <- rep(1,C)
  y2 <- as.vector(t(y))
  x2 <- t(kronecker(rep(1,C),x))
  dim(x2) <- c(r,n,C)
  x2 <- aperm(x2,c(1,3,2))
  dim(x2) <- c(r,n*C)
  x2 <- t(x2)

  X <- matrix(0,nrow=C*n,ncol=C*r)
  for ( i in 1:n)
    {
      for (j in 1:C)
        {
          X[((i-1)*C+j),(((j-1)*r+1):((j-1)*r+r))] <- x[i,]
        }
    }
  
  dx <- dim(X)
  
  if (is.null(K.prov)) K <- min(dx[1]-1,dx[2])
  else if ( K.prov > min(dx[1]-1,dx[2]) )
    {
      cat("number of dimension K.provd exceeds rank of X. \n","Provide a number that is less than or equal to ",min(dx[1]-1,dx[2]),"!\n"); break
    }
  else
    {
      K <- K.prov
    }

#  cat("Number of components is:", K,"!\n")

### initialzing matrices for PLS regression

  ## weight matrix for PLS regression
  W<- matrix(0,dx[2],K)

  ## score matrix
  Ta <- matrix(0,dx[1],K)

  ## loading matrix for X
  P <- matrix(0,dx[2],K)

  ## loading matrix for y
  Q <- numeric(K)


### intial values

  ## convergence
  converged <- F
  
  ## initial predictor matrix
  E <- X

  ## pseudo response
  ystar <- psi(y2,family)
  ystar0 <- ystar 

  ## Weight matrix for GLM
  V <- diag(abs(rnorm(C*n)),C*n)

  ## linear predictor, ie. eta = Xb
  eta.old <- rep(1,length(y2))
  eta <- eta.old+100

  Q <- rep(0,K)
  Q.old <- rep(100,K)
  K.old <- K
  min <- 1000
  
  ## regression coefficients
  if(!is.null(b.ini))
    {
      beta <- b.ini
    }
  else{
    beta <- rep(1,C*r)
  }

  beta.old <- beta/1000

  ## iteration index
  l <- 0
  l.min <- 1

  while ((max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)) > eps) & (l < lmax))
    {
  #    cat("Iter:",l,"\n")

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
    
      temp <- ginv(t(P)%*%W)
      if(any(is.na(temp)))
        {
          #cat("t(P)%*%W singular!\n")
          break
        }
      else
        {
          beta.old <- beta
          beta <- W%*%temp%*%Q
        }
    #   beta <- W%*%ginv(t(P)%*%W)%*%Q  

      ## Hat matrix
      H <- V%*%X%*%ginv(t(X)%*%V%*%X)%*%t(X)

      ## update eta (linear predictor)

      eta <- as.vector(t(x%*%matrix(beta,ncol=C,byrow=F)))
      p <-exp(x%*%matrix(beta,ncol=C ,byrow=F))
     
      den.p <- 1+as.vector(p%*%jC)
      
       if(any(is.infinite(den.p)))
        {
          #cat("Infinite denom for p!\n")
          break
        }
      p <- p2 <- p/den.p
      dim(p2) <- c(n,C)
      p2 <- as.vector(t(p2))

      ## update and rescaling weight matrix
      V <- matrix(0, C*n, C*n)
      for(j in seq(n)) {
        
        V[((j - 1) * C + 1):(j * C), ((j - 1) * C + 1):(j * C)] <- diag(p[j,],length(p[j,]))-matrix(p[j,],nrow=C,ncol=1)%*%matrix(p[j,],nrow=1,ncol=C)

      }
    
       ## diagnosis for divergence
      
    ## d_eta/d_mu matrix
      detadmu <- V
     	 for ( j in seq(n))
           {
             for (i in ((j - 1) * C + 1):(j * C))
               for (k in (i:(j*C)))
                 {
                   sum <- sum(p2[((j-1)*C+1):(j*C)])
                   if(i==k)
                     {
                     
                       detadmu[i,k] <- (1-sum+p2[k])/(p2[k]*(1-sum))
                         
                     }
                   else
                     {
                       detadmu[i,k] <- 1/(1-sum)
                       detadmu[k,i] <- detadmu[i,k]
                     }

                 }
           }
      
      Hw <- NULL
      for ( j in seq(n))
        {
          sum <- sum(diag(H)[((j - 1) * C + 1):(j * C)])
          for (i in ((j - 1) * C + 1):(j * C)) 
            {
              Hw[i] <- sum + diag(H)[i]
            }
        }
      
    ## update pseudo response
      if(any(is.na(detadmu)))
        {
          #cat("V singular\n")
          break
        }
     
      ystar <- eta + detadmu%*%(y2+diag(H)*br/2-(Hw*br/2+1)*p2)/(Hw*br/2+1)
      ystar0 <- ystar
      V <- V*(Hw*br/2+1)

      ## update predictor matrix
      E <- X
     
#      print(max(abs(beta-beta.old)/abs(beta.old)))
    }

 # print(l)
  if (max(abs(beta-beta.old)/(abs(beta.old)+denom.eps)) > eps)
    {
      
   #   cat("Convergence Not achieved and estimate from iteration ", l.min, " is used!\n")
      W <- W.min
      P <- P.min
      Q <- Q.min
    }
  else
    {
      converged <- T
    }
  ## final estimates
  beta <- W%*%ginv(t(P)%*%W)%*%Q

return( list(coefficients=matrix(beta,ncol=C ,byrow=F),     
              convergence = converged,
              niter = l,bias.reduction =br))
       
}
