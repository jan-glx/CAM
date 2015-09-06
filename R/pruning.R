pruning <-
function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) 
{
    p <- dim(G)[1]
    finalG <- sparseMatrix(i=integer(),j=integer(),dims=c(p,p))
    for(i in 1:p)
    {
        parents <- which(G[,i])
        lenpa <- length(parents)

        if(output)
        {
            cat("pruning variable:", i, "\n")
            cat("considered parents:", parents, "\n")
        }
        
        if(lenpa>0)
        {
            Xtmp <- cbind(X[,parents],X[,i])
            selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
            finalParents <- parents[selectedPar]
            finalG[finalParents,i] <- TRUE
        }
    }
    
    return(finalG)
}
