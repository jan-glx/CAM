#' @
#' @export
fastForward.cam <- function(object, noise)
{   
    output <- copy(noise) # copy unneccesarly as noise is not used elsewhere, but lets be friendly and have no sideeffects
    causalOrder <- dagToCausalOrder(object$causalDAG)
    for (k in causalOrder){
        set(output,i=NULL , j=k, value=predict(object$nodeModels[[k]], output) + noise[[k]])
    }
    return(output)
}

#' @export
colwise_resample <- function(X, seed_=NULL){
    if (!is.null(seed_)) {seed.bak <- .GlobalEnv$.Random.seed; set.seed(seed_)}
    ret <- setDT(lapply(X, function(x) x[sample(length(x), length(x), replace=TRUE)]))
    if (!is.null(seed_)) .GlobalEnv$.Random.seed <- seed.bak
    return(ret)
}

#' Bootstrap!
#' @param method either of "two-sided" or "one-sided"
#' @export
#'
#' @examples
#'
#' p=5
#' trueDAG <- matrix(FALSE,ncol=p, nrow=p)
#' trueDAG[matrix(c(1,2,3,4,
#'                  3,4,5,5),ncol=2)] <- TRUE
#' obj <- random_additive_polynomial_SEM(trueDAG, degree=2, seed_=1)
#' X <- simulate_additive_SEM(obj, n=100, seed_=3)
#' boot_res <- bootstrap.cam(X, matrix(c(2,1,1,2), nrow=2, byrow = TRUE),B=20) #two-sided null
#' boot_res$pvalue
bootstrap.cam <- function(X, fixedOrders, B=100, scoreName="SEMLINPOLY", parametric=FALSE, 
                          bootstrapH02 = FALSE, fast_double=FALSE, parsScore=list(), verbose=0,
                          method = NULL){
    X <- as.data.table(X)
    fixedOrdersH0 <- fixedOrders
    fixedOrdersOther <- matrix(nrow=0,ncol=2)
    if(!is.null(method) && method=="two-sided") {
        fixedOrdersH0 <- rbind(fixedOrdersH0,fixedOrdersH0[,c(2,1)])
    }else if(!is.null(method) && method=="one-sided") {
        fixedOrdersOther <- fixedOrdersH0[,c(2,1)]
    } else {
        warning("unknown method: \"",method,"\"-leaving fixedOrders as is.")
    }
    full_model <- CAM(X, orderFixationMethod="emulate_edge", fixedOrders=fixedOrdersOther, scoreName=scoreName, parsScore=parsScore)
    null_model <- CAM(X, orderFixationMethod="emulate_edge", fixedOrders=fixedOrdersH0, 
                      scoreName= scoreName,parsScore=parsScore)
    null_fit <- cam.fit(X, null_model$Adj, scoreName=scoreName,parsScore=parsScore)
    resid <- residuals(null_fit)
    log_likelihood_ratio <- full_model$Score-null_model$Score

    if (parametric){
      sds <- lapply(resid, sd)
    }

    if (fast_double){
        second_level_stats <- rep(NA,B)
    }

    stat_matrix <- matrix(NA,ncol=3, nrow=B)
    for (b in 1:B){
        if(verbose) print(b)
        if (bootstrapH02) {
            Xboot <- X[sample(nrow(X), nrow(X), replace=TRUE)]
            boot1_full_model <- CAM(Xboot, orderFixationMethod="emulate_edge", fixedOrders=fixedOrdersOther, scoreName=scoreName, parsScore=parsScore)
            boot1_null_model <- CAM(Xboot, orderFixationMethod="emulate_edge", fixedOrders=fixedOrdersH0, 
                              scoreName= scoreName,parsScore=parsScore)
            null_fit <- cam.fit(Xboot, boot1_null_model$Adj, scoreName=scoreName,parsScore=parsScore)
            resid <- residuals(null_fit)
        }
            
        if (!parametric) {
          resid_boot <- colwise_resample(resid)
        } else {
          resid_boot <- setDT(lapply(sds, function(x) rnorm(nrow(X), 0, x)))
        }
        Xboot <- fastForward.cam(null_fit, resid_boot)

        boot_full <-  CAM(Xboot, orderFixationMethod="emulate_edge", fixedOrders=fixedOrdersOther, scoreName= scoreName, parsScore=parsScore)
        boot_null <-  CAM(Xboot, orderFixationMethod="emulate_edge",
                          fixedOrders=fixedOrdersH0, scoreName= scoreName, parsScore=parsScore)

        stat_matrix[b,1] <- boot_full$Score
        stat_matrix[b,2] <- boot_null$Score

        if (fast_double){
            second_level_stats[b] <- bootstrap.cam(Xboot, fixedOrdersH0, B=1, scoreName=scoreName, 
                parametric=parametric, parsScore=parsScore)$bootstrap_samples[1,3]
        }
    }
    stat_matrix[,3] <- stat_matrix[,1]-stat_matrix[,2]
    pvalue <- mean(log_likelihood_ratio < stat_matrix[,3])+ mean(log_likelihood_ratio == stat_matrix[,3])/2 #mid pvalue 1. Franck, W. E. P-values For Discrete Test Statistics. Biometrical J. 28, 403–406 (1986).
    res <- list(log_likelihood_ratio=log_likelihood_ratio,
        pvalue=pvalue,
        bootstrap_samples = stat_matrix, 
        full_model=full_model$Score, 
        null_model=null_model$Score,
        parametric=parametric, 
        bootstrapH02 = bootstrapH02, 
        fast_double = fast_double,
        method = method
        )

    if (fast_double){
        quant <- quantile(second_level_stats, 1-pvalue)
        double_pvalue <- mean(quant < stat_matrix[,3])
        res$double_pvalue <- double_pvalue
    }
    return(res)
}