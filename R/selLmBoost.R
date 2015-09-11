#' auxiliary file for CAM: selection based on lm (mboost)
#' 
#' auxiliary file for CAM: selection based on lm (mboost)
#' 
#' 
#' @export selLmBoost
selLmBoost <-
function(X,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE,k)
{
    return(selXXBoost(train_LMboost, X, pars, output, k))
}
