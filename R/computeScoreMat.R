#' auxiliary file for CAM: Computes the initial score matrix.
#' 
#' auxiliary file to CAM. Computes the initial score matrix.
#' 
#' 
#' @param X nxp matrix of training inputs (n data points, p dimensions)
#' @param scoreName specifies the model type which is used to compute the
#' score. Default is "SEMGAM" which assumes a generalized additive model class.
#' Other options include "SEMLIN" which fits a linear model.
#' @param numParents indicates how many parents we consider. If numParents = 1
#' (default), then the score matrix is of dimension (p-1) x p. If numParents =
#' 2, then the score matrix is of dimension (p-1)(p-2) x p and so on
#' @param output boolean indicating whether information about the progress is
#' written to the console.
#' @param numCores specifies the number of cores that can be used for
#' computation.
#' @param selMat indicating the possible parent relationships.
#' @param parsScore additional parameters can be supported to the score
#' function.
#' @param intervMat the matrix intervMat has the same dimension as X. entry
#' (i,j) == TRUE indicates that in experiment i, variable j has been intervened
#' on.
#' @param intervData boolean that indicates whether we use interventional data.
#' @return A list with elements \item{scoreMat}{The score matrix. scoreMat[i,j]
#' contains the gain in score if we consider i being a parent of j }
#' \item{rowParents}{Contains the row names of the score matrix. Only relevant
#' if numParents > 1.} \item{scoreEmtpyNodes}{Vector containing the scores of
#' each node in the empty graph without any edges.}
#' @note This is an auxiliary file for CAM.
#' @author J. Peters (jonas.peters@@tuebingen.mpg.de) and J. Ernest
#' (ernest@@stat.math.ethz.ch)
#' @seealso \code{\link[CAM]{CAM}}
#' @references P. B\"uhlmann, J. Peters, J. Ernest: CAM: Causal Additive
#' Models, high-dimensional Order Search and Penalized Regression Annals of
#' Statistics 42:2526-2556, 2014.
#' @export computeScoreMat
#' @import data.table
computeScoreMat <-
function(X, scoreName, numParents, output, numCores, selMat, parsScore, intervMat, intervData)
{
    
    # numParents indicates how many parents we consider. If numParents = 1 (default), then the 
    # score matrix is of dimension (p-1) x p. If numParents = 2, then the  
    # score matrix is of dimension (p-1)(p-2) x p and so on...
    #
    # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
    # it should therefore be positive.
    
    p <- dim(X)[2]
    n <- dim(X)[1]
    rowParents <- t(combn(p,numParents))
    
    tt <- expand.grid(1:dim(rowParents)[1], 1:p)
    allNode2 <- tt[,2]
    allI <- tt[,1]
    argList <- list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, output = output,
                    parsScore = parsScore, intervMat = intervMat, intervData = intervData)
    if(numCores == 1)
    {
        scoreMat <- mapply(computeScoreMatParallel, MoreArgs = argList, node2 = allNode2, i = allI)
    } else
    {
        scoreMat <- parallel::mcmapply(computeScoreMatParallel, MoreArgs = argList, node2 = allNode2, i = allI, mc.cores = numCores)
    }
    
    scoreMat <- matrix(scoreMat,dim(rowParents)[1],p)
    # initScore[i] equals the variance of variable j. 
    initScore <- rep(NA,p)
    for(i in 1:p)
    {
        if(intervData)
        {
            X2 <- X[!intervMat[,i],]
        } else
        {
            X2 <- X
        }
        vartmp <- var(X2[,i])
        initScore[i] <- -log(vartmp)
        # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
        scoreMat[,i] <- scoreMat[,i] - initScore[i]
    }
    return(list(scoreMat = scoreMat, rowParents = rowParents, scoreEmptyNodes = initScore))
}
