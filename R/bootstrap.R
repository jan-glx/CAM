#' @export
fastForward.cam <- function(object, noise)
  {   output <- as.data.table(noise)
  p <- object$p
  causalOrder <- CAM::dagToCausalOrder(as.matrix(object$causalDAG))
  for (k in causalOrder){
    output[,k] <- predict(object$nodeModels[[k]], output) + noise[,k]
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
#'                 4,5),ncol=2,byrow=T)] <- TRUE
#' obj <- random_additive_polynomial_SEM(trueDAG, degree=2, seed_=1)
#' X <- simulate_additive_SEM(obj, n=100, seed_=3)
#' boot_res <- bootstrap.cam(X, matrix(c(2,1,1,2), nrow=2, byrow = TRUE),B=20) #two-sided null
#' boot_res$pvalue

bootstrap.cam <- function(X, fixedOrders, B=100, scoreName="SEMLINPOLY"){
    full_model <- CAM(X, pruning=FALSE,orderFixationMethod="emulate_edge",scoreName=scoreName)
    null_model <- CAM(X, pruning=FALSE,orderFixationMethod="emulate_edge",
                  fixedOrders=fixedOrders, scoreName= scoreName)
    null_fit <- cam.fit(X, null_model$Adj)
    resid <- CAM:::residuals.cam(null_fit)
    stat_matrix <- matrix(NA,ncol=3, nrow=B)
    for (b in 1:B){
        print(b)
        resid_boot <- colwise_resample(resid)
        Xboot <- fastForward.cam(null_fit, resid_boot)
        boot_full <-  CAM(Xboot, pruning=FALSE,orderFixationMethod="emulate_edge",scoreName= scoreName)$Score
        boot_null <-  CAM(Xboot, pruning=FALSE,orderFixationMethod="emulate_edge",
                          fixedOrders=fixedOrders, scoreName= scoreName)$Score
        stat_matrix[b,1] <- boot_full
        stat_matrix[b,2] <- boot_null
    }
    stat_matrix[,3] <- stat_matrix[,1]-stat_matrix[,2]
    log_likelihood_ratio <- full_model$Score-null_model$Score
    pvalue <- mean(log_likelihood_ratio < stat_matrix[,3])
    return(list(log_likelihood_ratio=log_likelihood_ratio,
        pvalue=pvalue,
        bootstrap_samples = stat_matrix, 
        full_model=full_model$Score, 
        null_model=null_model$Score))
}