selXXBoost <-
  function(train_function, X, pars, output, k)
  {
    if(output)
    {
      cat("Performing variable selection for variable", k, ": \n")
    }
    p <- ncol(X)
    selVec <- rep(FALSE, p)
    if(p >= 2)
    {
      modfitGam <- train_function(X[,-k], X[,k], pars)
      dt <- data.table(variable = modfitGam$model$xselect())
      dt <- dt[,.N,by=variable][,.(variable, N=N/nrow(dt))][N>0,]
      # ^^^ O(n*log(p)), O(n) possible as we know possible values
      if(output)
      {
        cat("The following variables \n")
        show(dt$variable)
        cat("... have been selected in the following fractions: \n")
        show(dt$N)
      }
      dt <- dt[N>=pars$atLeastThatMuchSelected]
      dt <- dt[order(N)[0:min(pars$atMostThatManyNeighbors,nrow(dt))]]
      if(output)
      {
        cat("We finally choose as possible parents: \n")
        show(dt$variable)
        cat("\n")
      }
      tmp <- rep(FALSE,p-1)
      tmp[dt$variable] <- TRUE
      selVec[-k] <- tmp
    }
    return(selVec)
  }
