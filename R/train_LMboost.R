train_LMboost <-
function(X,y,pars = list()) #
{
    # centralize X
    #X <- X - matrix( rep( colSums(X)/dim(X)[1],dim(X)[1]), dim(X)[1], dim(X)[2], byrow = TRUE) 
    y <- y - rep( mean(y), length(y))
    
    yy <- as.vector(y)
    op <- options(warn=-1)
    on.exit(options(op))
    gb <- glmboost(X,yy, center = TRUE)
    
    result <- list()
    result$Yfit <- gb$fitted()
    result$residuals <- gb$resid()
    result$model <- gb
    return(result)
}
