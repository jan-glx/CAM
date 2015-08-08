selLmBoost <-
function(X,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE,k)
{
    return(selXXBoost(train_LMboost, X, pars, output, k))
}
