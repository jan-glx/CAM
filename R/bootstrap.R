#' @export
fastForward.cam <- function(object, noise)
{   
    noise <- as.data.table(noise)
    output <- noise
    p <- object$p
    causalOrder <- dagToCausalOrder(as.matrix(object$causalDAG))
    for (k in causalOrder){
        if (any(object$causalDAG[,k])){
            output[,k] <- predict(object$nodeModels[[k]], output) + noise[,k,with=F][[1]]
        }
    }
    as.matrix(output)
}

#' @export
colwise_resample <- function(X, seed_=NULL){
    if (!is.null(seed_)) set.seed(seed_)
    apply(X, 2, function(x) x[sample(1:length(x), length(x), replace=TRUE)])
}

#' Bootstrap!
#' @export
#'
#' @examples
#'
#' p=5
#' trueDAG <- matrix(FALSE,ncol=p, nrow=p)
#' trueDAG[matrix(c(1,3,
#'                  2,4,
#'                  3,5,
#'                  4,5),ncol=2,byrow=TRUE)] <- TRUE
#' obj <- random_additive_polynomial_SEM(trueDAG, degree=2, seed_=1)
#' X <- simulate_additive_SEM(obj, n=100, seed_=3)
#' boot_res <- bootstrap.cam(X, matrix(c(2,1,1,2), nrow=2, byrow = TRUE),B=20) #two-sided null
#' boot_res$pvalue

bootstrap.cam <- function(X, fixedOrders, B=100, scoreName="SEMLINPOLY", parametric=FALSE, fast_double=FALSE, parsScore=list()){
    full_model <- CAM(X, scoreName=scoreName, parsScore=parsScore)
    null_model <- CAM(X, orderFixationMethod="emulate_edge", fixedOrders=fixedOrders, 
                      scoreName= scoreName,parsScore=parsScore)
    null_fit <- cam.fit(X, null_model$Adj, scoreName=scoreName,parsScore=parsScore)
    resid <- residuals.cam(null_fit)

    if (parametric){
      sds <- apply(resid, 2, sd)
    }

    if (fast_double){
        second_level_stats <- rep(NA,B)
    }

    stat_matrix <- matrix(NA,ncol=3, nrow=B)
    for (b in 1:B){
        print(b)
        if (!parametric){
          resid_boot <- colwise_resample(resid)
        } else {
          resid_boot <- sapply(sds, function(x) rnorm(nrow(X), 0, x))
        }
        Xboot <- fastForward.cam(null_fit, resid_boot)

        boot_full <-  CAM(Xboot, scoreName= scoreName, parsScore=parsScore)
        boot_null <-  CAM(Xboot, orderFixationMethod="emulate_edge",
                          fixedOrders=fixedOrders, scoreName= scoreName, parsScore=parsScore)

        stat_matrix[b,1] <- boot_full$Score
        stat_matrix[b,2] <- boot_null$Score

        if (fast_double){
            second_level_stats[b] <- bootstrap.cam(Xboot, fixedOrders, B=1, scoreName=scoreName, 
                parametric=parametric, parsScore=parsScore)$bootstrap_samples[1,3]
        }
    }
    stat_matrix[,3] <- stat_matrix[,1]-stat_matrix[,2]
    log_likelihood_ratio <- full_model$Score-null_model$Score
    pvalue <- mean(log_likelihood_ratio < stat_matrix[,3])
    res <- list(log_likelihood_ratio=log_likelihood_ratio,
        pvalue=pvalue,
        bootstrap_samples = stat_matrix, 
        full_model=full_model$Score, 
        null_model=null_model$Score)

    if (fast_double){
        quant <- quantile(second_level_stats, 1-pvalue)
        double_pvalue <- mean(quant < stat_matrix[,3])
        res$double_pvalue <- double_pvalue
    }
    return(res)
}