glpls1a.train.test.error <- function(train.X, train.y, test.X, test.y,
                                     K.prov=NULL, eps=1e-3, lmax=100,
                                     family="binomial", link="logit", br=T)
{
  ## calculate test set error by out-of-sample counting
  ## the input predictor matrices are p by n

  train.X <- t(train.X)
  test.X <- t(test.X)

  if(!is.null(rownames(train.X)) & !is.null(rownames(test.X)))
     {
        test.X <- test.X[match(row.names(train.X),row.names(test.X)),]
     }

  train.mean <- rowMeans(train.X)
  train.sd <- apply(train.X,1,sd)
  train.X <- scale(t(train.X))

  test.X <- t((test.X-train.mean)/train.sd)

  train.fit <- glpls1a(train.X, train.y, K.prov=K.prov, eps=eps,
                      lmax=lmax, family=family, link=link, br=br)
  test.y.predict <- glpls1a.predict(test.X, train.fit$coefficients,
                                       family=family, link=link)
  error.test <- 1-sum((test.y.predict>0.5)==test.y)/length(test.y)
  error.obs <- seq(1:length(test.y))[(test.y.predict>0.5)!=test.y]
  return(list(error=error.test, error.obs = error.obs,
              predict.test=test.y.predict))

}

glpls1a.cv.error <- function(train.X, train.y, K.prov=NULL, eps=1e-3,
lmax=100, family="binomial", link="logit", br=T)
  {
    ## calculate error rate by leaving-one-out CV of training set
    ## the input train.X is p by n

    train.y.predict <- NULL
    for( i in 1:nrow(train.X))
      {
    #    print(i)
        train.fit <- glpls1a(train.X[-i,], train.y[-i], K.prov=K.prov,
        eps=eps,lmax=lmax,family=family,link=link,br=br)
        train.y.predict[i] <-
        glpls1a.predict(matrix(train.X[i,], nrow=1), train.fit$coefficients)
      }

    error.CV <- 1-sum((train.y.predict>0.5)==train.y)/length(train.y)
    error.obs <- seq(1:length(train.y))[(train.y.predict>0.5)!=train.y]
    return(list(error=error.CV,error.obs=error.obs))
  }

glpls1a.error <- function(train.X, train.y, K.prov=NULL, eps=1e-3,
lmax=100, family="binomial", link="logit", br=T)
  {
    ## calculate error rate by leaving-one-out CV of training set
    ## the input train.X is p by n

    train.y.predict <- NULL
    train.fit <- glpls1a(train.X, train.y, K.prov=K.prov, eps=eps,
    lmax=lmax, family=family, link=link, br=br)
    train.y.predict <- glpls1a.predict(train.X, train.fit$coefficients)
    error <- 1-sum((train.y.predict>0.5)==train.y)/length(train.y)
    error.obs <- seq(1:length(train.y))[(train.y.predict>0.5)!=train.y]
    return(list(error=error,error.obs=error.obs))
  }


glpls1a.mlogit.cv.error <- function(train.X, train.y, K.prov=NULL,
eps=1e-3, lmax=100, mlogit=T, br=T)
  {
    ## calculate CV error rate for multinomial logit model
    ## calculate error rate by leaving-one-out CV of training set
    ## the input train.X is p by n
    ## the input train.y is the multinomial outcome, within level 1 being the most frequent

    family <- "binomial"
    link <- "logit"

    n <- nrow(train.X)

    C <- max(train.y)
    train.y.predict <- matrix(0,nrow=n,ncol=C-1)

    if(is.null(K.prov))
      {
        dx <- dim(train.X)
        K.prov <- min(dx[1],dx[2]-1)
      }


    for( i in 1:n)
      {
    #   print(i)
       if(mlogit)
         {
           train.fit <- glpls1a.mlogit(train.X[-i,], train.y[-i],
           K.prov=K.prov, eps=eps, lmax=lmax, br=br)
         }
       else
         {
            train.fit <- glpls1a.logit.all(train.X[-i,], train.y[-i],
            K.prov=K.prov, eps=eps, lmax=lmax, br=br)
          }

       train.y.predict[i,] <-
       glpls1a.mlogit.predict(matrix(train.X[i,],nrow=1), train.fit$coef)
      }

    temp <- t(apply(train.y.predict,1,function(x) x==max(x,1-sum(x))))
    temp <- apply(cbind(apply(temp,1,function(x) !any(x)),temp), 1,
                                        function(x) (1:C)[x])
    error.CV <- 1-sum(temp==train.y)/length(train.y)
    error.obs <- seq(1:length(train.y))[temp != train.y]
    return(list(error=error.CV,error.obs=error.obs))
  }
