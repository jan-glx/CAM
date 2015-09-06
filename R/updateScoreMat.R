updateScoreMat <-
function(scoreMat, X, scoreName, i, j, scoreNodes, Adj, output, numCores, maxNumParents, parsScore, intervMat, intervData)
    # new edge: from i to j
{
    p <- dim(X)[2]
    existingParOfJ <- which(Adj[,j])
    notAllowedParOfJ <- which(scoreMat[,j] == -Inf)
    # if there is something left that we need to update
    if(length(notAllowedParOfJ) < p)
    {
        # update column for j
        rowParents <- matrix(c(existingParOfJ,NA), p, length(existingParOfJ)+1, byrow = TRUE)
        rowParents[,length(existingParOfJ)+1] <- 1:p
        toUpdate <- setdiff(1:p,notAllowedParOfJ)
        if(length(existingParOfJ)< maxNumParents)
        {
            argList <- list(rowParents = rowParents, selMat = matrix(TRUE,p,p), 
                            scoreName = scoreName, X = X, output = output, node2 = j, 
                            parsScore = parsScore, intervMat = intervMat, intervData = intervData)
            if(numCores == 1)
            {
                scoreUpdate <- mapply(computeScoreMatParallel,MoreArgs = argList, i = toUpdate)
            } else
            {
                scoreUpdate <- mcmapply(computeScoreMatParallel,MoreArgs = argList, i = toUpdate, mc.cores = numCores)
            }
        } else
        {
            scoreUpdate <- -Inf
        }
        scoreMat[toUpdate,j] <- scoreUpdate - scoreNodes[j]
    }    
    return(scoreMat)
}
