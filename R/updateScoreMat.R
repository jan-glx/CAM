#' auxiliary file for CAM: Updates the score matrix after having added edge i
#' -> j to the graph.
#' 
#' auxiliary file for CAM. Updates the score matrix after having added edge i
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
#' @param scoreNodes vector containing current scores of each node
#' @param Adj adjacency matrix of the graph
#' @param output boolean indicating whether information about the progress is
#' written to the console.
#' @param numCores specifies the number of cores that can be used for
#' computation.
#' @param maxNumParents specifies the maximal number of parents that are
#' allowed in the model.
#' @param parsScore additional parameters can be supported to the score
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
                scoreUpdate <- parallel::mcmapply(computeScoreMatParallel,MoreArgs = argList, i = toUpdate, mc.cores = numCores)
            }
        } else
        {
            scoreUpdate <- -Inf
        }
        scoreMat[toUpdate,j] <- scoreUpdate - scoreNodes[j]
    }    
    return(scoreMat)
}
