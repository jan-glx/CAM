train_additive_polynomial <-
function(X,y,pars = list(degree = 3))
{
    if(!("degree" %in% names(pars) ))
    { 
        pars$degree = 3
    }
    p <- dim(as.matrix(X))

    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
        coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
        for(i in 2:p[2])
        {
            labs<-paste(labs,"poly(var",i,",degree = ",pars$degree,") + ",sep="")
        }
    }
    labs<-paste(labs,"poly(var",p[2]+1,",degree = ",pars$degree,")",sep="")
    mod_additive_polynomial <- lm(formula=formula(labs), data=dat)
    
    result <- list()
    result$Yfit <- as.matrix(mod_additive_polynomial$fitted.values)
    result$residuals <- as.matrix(mod_additive_polynomial$residuals)
    result$model <- mod_additive_polynomial 
    result$df <- mod_additive_polynomial$df.residual     
    #result$edf <- mod_additive_polynomial$edf     
    #result$edf1 <- mod_additive_polynomial$edf1     
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
}
