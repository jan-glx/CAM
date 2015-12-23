

incEdge <- function(X, nodeModelName, nodeModelPars, scoreFunction, maxNumParents, fixedOrders, 
                    orderFixationMethod, selMat, intervData, intervMat, numCores= 1, verbose = FALSE){
    
    if (is.list(X)||is.data.frame(X)||is.data.table(X)) setDT(X)
    else X <- as.data.table(X)
    
    # we record how the score develops 
    scoreVec <- c()
    # and which edges are added
    edgeList <- c()
    
    p <- ncol(X)
    
    # this counter is only used if verbose = TRUE
    counterUpdate <- 0
    # We need the pathMatrix (entry (i,j) being one means that there is a directed path from i to j) in order to keep track of possible cycles.
    pathMatrix <- matrix(FALSE,p,p)
    diag(pathMatrix) <- TRUE
    if (orderFixationMethod == "emulate_edge"){
        fixedOrdersMatrix <- matrix(FALSE, p, p)
        fixedOrdersMatrix[fixedOrders] <- TRUE
        # prevent addition of edges from DE(j) to {AC(i),i} and {DE(j),j} to {AC(i)} where (i,j) is a row of fixedOrders
        pathMatrix[getPathMatrix(fixedOrdersMatrix)] <- TRUE
    }
    scoreMat <- matrix(NA_real_, p, p)
    scoreMat[!selMat] <- -Inf
    scoreMat[t(pathMatrix)] <- -Inf
    Adj <- matrix(FALSE, p, p)
    
    # initialize score matrix
    intercept_only_model <- cam.fit(X, causalDAG=Adj, nodeModelName, nodeModelPars, numCores=numCores)
    scoreNodes <- sapply(intercept_only_model$nodeModels, scoreFunction)
    for (j in 1:p) {
        scoreMat <- updateScoreMat(scoreMat=scoreMat, X=X, nodeModelName=nodeModelName, 
                                   nodeModelPars=nodeModelPars, scoreFunction=scoreFunction, i=NULL, 
                                   j=j, scoreNodes=scoreNodes, Adj=Adj, maxNumParents=maxNumParents, 
                                   intervMat=intervMat, intervData=intervData, verbose=verbose, 
                                   numCores=numCores)
    }

    fixedOrdersAdded <- 0L
    # Greedily adding edges
    while(sum(scoreMat!=-Inf) > 0)
    {
        if (orderFixationMethod == "force_edge" && fixedOrdersAdded < nrow(fixedOrders)){
            fixedOrdersAdded = fixedOrdersAdded + 1L
            ix_max <- fixedOrders[fixedOrdersAdded, , drop=F]
        } else {
            ix_max <- arrayInd(which.max(scoreMat), dim(scoreMat))
        }
        
        Adj[ix_max] <- TRUE
        scoreNodes[ix_max[2]] <- scoreNodes[ix_max[2]] + scoreMat[ix_max]
        if(verbose)
        {
            cat("\n Included edge (from, to) ", ix_max, "\n")
        }
        
        # Do not include the same edge twice.
        scoreMat[ix_max] <- -Inf
        
        # Avoid cycles
        DescOfAndNewChild <- which(pathMatrix[ix_max[2],])
        AncOfAndNewParent <- which(pathMatrix[,ix_max[1]])
        pathMatrix[AncOfAndNewParent,DescOfAndNewChild] <- TRUE
        scoreMat[DescOfAndNewChild,AncOfAndNewParent] <- -Inf 
        
        # Record the score of the current graph
        scoreVec <- c(scoreVec, sum(scoreNodes))
        # Record which edge has been added
        edgeList <- rbind(edgeList, ix_max, deparse.level=0)
        
        # Update column j of scoreMatrix
        scoreMat <- updateScoreMat(scoreMat=scoreMat, X=X, nodeModelName=nodeModelName, 
                                   nodeModelPars=nodeModelPars, scoreFunction=scoreFunction, 
                                   i=ix_max[1], j=ix_max[2], scoreNodes=scoreNodes, Adj=Adj, 
                                   maxNumParents=maxNumParents, intervMat=intervMat, 
                                   intervData=intervData, verbose=verbose, numCores=numCores)
        counterUpdate <- counterUpdate + 1
    }
    return(list(Adj = Adj, scoreNodes=scoreNodes, scoreVec = scoreVec, edgeList = edgeList))
}

#' Computes an entry of the score matrix that is the gain in scoreFunctions value upon addition of an edge from i to j
#' 
#' @note This is an auxiliary file for the function CAM in package CAM.
#' @author Jonas Peters <jonas.peters@@tuebingen.mpg.de> and Jan Ernest
#' <ernest@@stat.math.ethz.ch>
#' @seealso \code{\link[CAM]{CAM}}
#' @references P. B\"uhlmann, J. Peters, J. Ernest: CAM: Causal Additive
#' Models, high-dimensional Order Search and Penalized Regression Annals of
#' Statistics 42:2526-2556, 2014.
#' @export
computeScoreMatEntry <- function(existingParOfJ, nodeModelName, nodeModelPars, scoreFunction, X, 
                                    verbose, i, j, oldScore) {
    newParentsOfJ <- c(existingParOfJ, i)
    if(verbose) cat("\r compute score entry for regressing",j,"on",newParentsOfJ,"                  \r")
    if(j %in% newParentsOfJ) stop(paste0("j(",j,") can not be a parent of itself. Ensure that i!=j and ! j %in% existingParOfJ!"))
    
    nodeModel <-fitNode(X, j, parents_of_j = newParentsOfJ, method = nodeModelName, pars = nodeModelPars) 
    score <- scoreFunction(nodeModel)
    return(score-oldScore)
}
#' auxiliary file for CAM: Updates the score matrix after having added edge i
#' -> j to the graph.
#' 
#' 
#' @param scoreMat the current score matrix that has to be updated
#' @param X nxp matrix of training inputs (n data points, p dimensions)
#' @param scoreName specifies the model type which is used to compute the
#' score. Default is "SEMGAM" which assumes a generalized additive model class.
#' Other options include "SEMLIN" which fits a linear model.
#' @param i starting point of the edge i->j that has been added
#' @param j end point of the edge i->j that has been added
#' @param nodeModelName vector containing current scores of each node
#' @param Adj adjacency matrix of the graph
#' @param verbose boolean indicating whether information about the progress is
#' written to the console.
#' @param numCores specifies the number of cores that can be used for
#' computation.
#' @param maxNumParents specifies the maximal number of parents that are
#' allowed in the model.
#' @param nodeModelPars additional parameters can be supported to the score
#' function.
#' @param intervMat the matrix intervMat has the same dimension as X. entry
#' (i,j) == TRUE indicates that in experiment i, variable j has been intervened
#' on.
#' @param intervData boolean that indicates whether we use interventional data.
#' @return \item{scoreMat}{the updated score matrix.}
#' @note This is an auxiliary file for CAM.
#' @author J. Peters (jonas.peters@@tuebingen.mpg.de) and J. Ernest
#' (ernest@@stat.math.ethz.ch)
#' @seealso \code{\link[CAM]{CAM}}
#' @references P. B\"uhlmann, J. Peters, J. Ernest: CAM: Causal Additive
#' Models, high-dimensional Order Search and Penalized Regression Annals of
#' Statistics 42:2526-2556, 2014.
#' @export updateScoreMat
updateScoreMat <- function(scoreMat, X, nodeModelName, nodeModelPars, scoreFunction, i, j, 
                           scoreNodes, Adj, maxNumParents, intervMat, intervData, verbose = FALSE, 
                           numCores = 0) {
    # new edge: from i to j
    p <- dim(X)[2]
    existingParOfJ <- which(Adj[,j])
    if(length(existingParOfJ) >= maxNumParents) scoreMat[setdiff(1:p,existingParOfJ),j] <- -Inf
    notAllowedParOfJ <- which(scoreMat[,j] == -Inf)
    toUpdate <- setdiff(1:p,notAllowedParOfJ)
    if(length(toUpdate) > 0)
    {
        X2 <- if(intervData) X[!intervMat[,j],] else X
        argList <- list(existingParOfJ = existingParOfJ, nodeModelName = nodeModelName, 
                        X = X2, verbose = verbose, j = j, nodeModelPars = nodeModelPars, 
                        scoreFunction = scoreFunction, oldScore = scoreNodes[j])
        if(numCores == 1) {
            scoreUpdate <-             mapply(computeScoreMatEntry,MoreArgs = argList, i = toUpdate)
        } else {
            scoreUpdate <- parallel::mcmapply(computeScoreMatEntry,MoreArgs = argList, i = toUpdate, mc.cores = numCores)
        } 
        scoreMat[toUpdate,j] <- scoreUpdate 
    }
    return(scoreMat)
}
