glpls1a.predict <- function(X,beta,family="binomial",link="logit")
{
  if(all(X[,1] == rep(1,nrow(X))))
   {
     eta <- X%*%beta
   }else{
     eta <- cbind(rep(1,nrow(X)),X)%*%beta
  }
  return(h(eta,family,link))
}

glpls1a.mlogit.predict <- function(X,beta)
  {
    ## prediction for multinomial logit model
    ## beta is the p by J-1 matrix where J is the levels of categories of the outcomoe
    
    if(all(X[,1] == rep(1,nrow(X))))
      {
        eta <- X%*%beta
      }else
      {
        eta <- cbind(rep(1,nrow(X)),X)%*% beta
      }

    p <- exp(eta)
    p.denom <- apply(p,1,function(x) 1+sum(x))
      
    return(p/p.denom)
  }
