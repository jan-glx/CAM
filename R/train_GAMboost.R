train_GAMboost <-
function(X,y,pars = list()) #
{
    dat <- as.data.frame(X)
    bl <- lapply(dat, mboost::bbs)
    gb <- mboost::mboost_fit(bl, y)
    
    result <- list()
    result$Yfit <- gb$fitted()
    result$residuals <- gb$resid()
    result$model <- gb
    return(result)
}
