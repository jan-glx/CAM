selGam <-
function(X,pars = list(cutOffPVal = 0.001, numBasisFcts = 10),output = FALSE,k)
{
    if(output)
    {
      cat("Performing variable selection for variable", k, ": \n")
    }
    p <- ncol(X)
    selVec <- rep(FALSE, p)
    if(p >= 2)
    {
        mod_gam <- train_gam(X[,-k],as.matrix(X[,k]),pars)
        pValVec <- summary.gam(mod_gam$model)$s.pv
        if(output)
        {
            cat("vector of p-values:", pValVec, "\n")
        }
        stopifnot(length(pValVec) == length(selVec[-k]))
        selVec[-k] <- (pValVec < pars$cutOffPVal)
    }
    return(selVec)
}
