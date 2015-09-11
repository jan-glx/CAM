#' auxiliary file for CAM: selection based on gam boost (mboost)
#' 
#' auxiliary file for CAM: selection based on gam boost (mboost)
#' 
#' 
#' @export selGamBoost
selGamBoost <-
function(X,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE,k)
{
  return(selXXBoost(train_GAMboost, X, pars, output, k))
}
