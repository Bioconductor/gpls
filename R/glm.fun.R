"h" <-
  function(x,family="normal",link="identity") { 
# g-inverse, i.e. inverse funtion of link function
    if (family=="normal"){
      if (link=="identity") x
      else cat("link function not recognized for ",family, " choose from identity!\n")
    }
    else if (family=="binomial") {
      if (link=="logit")
        {
          ifelse(is.infinite(exp(x)),1,exp(x)/(1+exp(x)))
        
        }
      else if (link=="probit") pnorm(x)
      else if (link=="cloglog") 1-exp(-exp(x))
      else cat("link function not recognized for ",family, " choose from logit,probit and cloglog!\n")
    }
    else if (family=="poisson") {
      if (link=="log") exp(x)
      else cat("link function not recognized for ", family, " chosse from log!\n" )
    }
    else if (family=="gamma"){
      if (link=="reciprocal") 1/x
      else cat("link function not recognized for ",family, " choose from reciprocal!\n")
    }
  }


"hp" <-
function(x,family="normal",link="identity") { 
# first derivative of h wrt eta, the linear predictor

if (family=="normal"){
if (link=="identity") rep(1,length(x))
else cat("link function not recognized for ",family, " choose from identity!\n")}
else if (family=="binomial") {
if (link=="logit") exp(x)/(1+exp(x))^2
else if (link=="probit") dnorm(x)
else if (link=="cloglog") exp(-exp(x))*exp(x)
else cat("link function not recognized for ",family, " choose from logit,probit and cloglog!\n")
}
else if (family=="poisson") {
if (link=="log") exp(x)
else cat("link function not recognized for ", family, " chosse from log!\n" )
}
}

"g" <-
function(x,family="normal",link="identity") {
# mu-link function for generalized linear model
if (family=="normal"){
if (link=="identity") x
else cat("link function not recognized for ",family, " choose from identity!\n")}
if (family=="binomial") {
 if (link=="logit") log(x/(1-x))
 else if (link=="probit") qnorm(x)
 else if (link=="cloglog") log(-log(1-x))
 else cat("link function not recognized for ",family, " choose from logit,probitand cloglog!\n")
}
else if (family=="poisson") {
 if (link=="log") log(x)
 else cat("link function not recognized for ", family, " choose from log!\n" )
}
}

"bpp" <-
function(x,family="normal",link="identity") { 
# variance function
if (family=="normal"){
if (link=="identity") rep(1,length(x))
else cat("link function not recognized for ",family, " choose from identity!\n")}
else if (family=="binomial") {
if (link=="logit") exp(x)/(1+exp(x))^2
else if (link=="probit") pnorm(x)*(1-pnorm(x))
else if (link=="cloglog") exp(-exp(x))*(1-exp(-exp(x)))
else cat("link function not recognized for ",family, " choose from logit,probit and cloglog!\n")
}
else if (family=="poisson") {
if (link=="log") exp(x)
else cat("link function not recognized for ", family, " chosse from log!\n" )
}
}

"psi" <-
function(x,family="normal") {
# initial value for dependent variable
if (family=="binomial") return((x+0.5)/2)
else if (family=="poisson") return(log(x+0.5))
else return(x)}

"probcal" <- function(y,eta)
{ # return P(Y=y), now only works for logit link for binomial family
  p <- 0
  for ( i in 1:length(y))
  p[i] <- ifelse(y[i]==1, exp(eta[i])/(1+exp(eta[i])) ,1/(1+exp(eta[i])))
  p
}

