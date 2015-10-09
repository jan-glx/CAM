#' @export
random_additive_polynomial_SEM <- function(trueDAG, degree=3, noise_mean = 1, noise_variance = 0.5,
                                           intercept_variance = 1, rescale = TRUE, seed_=NULL) {
  if (!is.null(seed_)) {seed.bak <- .GlobalEnv$.Random.seed; set.seed(seed_)}
  p <- ncol(trueDAG)
  rand_poly <- function(.){
    coef <- runif(degree+1, -1, 1)
    f<- function(x){ t(sapply(x, `^`, (0:degree))) %*% coef}
    attr(f, "coef") <- coef
    return(f)
  }

  f_jk <- matrix(list(),p,p)
  f_jk[trueDAG] <- lapply(seq_len(sum(trueDAG)), FUN=rand_poly)
  mu <- rnorm(p)*intercept_variance
  e_var <- rnorm(p,noise_mean,noise_variance)
  result <- list(trueDAG = trueDAG, f_jk=  f_jk, mu = mu, p = p, e_var = e_var, scale = rep(1, p))
  if(rescale) result <- rescale_sem_object(result)
  if (!is.null(seed_)) .GlobalEnv$.Random.seed <- seed.bak
  return(result)
}

#' @export
rescale_sem_object <- function(sem_object, rsf = 3, n=500, seed_= NULL) {
    return(simulate_additive_SEM(sem_object, n, seed_= seed_, .rescale = TRUE, .rsf = rsf))
}



#' @export
simulate_additive_SEM <- function(sem_object, n = 500, seed_ = NULL, .rescale = FALSE, .rsf = 1) {
    if (!is.null(seed_)) {seed.bak <- .GlobalEnv$.Random.seed; set.seed(seed_)}
    
    p <- sem_object$p
    trueDAG <- sem_object$trueDAG
    f_jk <- sem_object$f_jk
    
    hsplit <- function(x) {
        x = as.matrix(x)
        split(x, col(x))
    }
    
    X <- matrix(rnorm(n * p), ncol = p) # noise
    for (k in order(dagToCausalOrder(sem_object$trueDAG))) {
        tmp  <- mapply(do.call,
                               f_jk[trueDAG[,k],k],
                               lapply(hsplit(X[,trueDAG[,k]]),list))
        if (.rescale & !(is.list(tmp))) {
            tmp2 <-  sem_object$mu[k] + 0 + 
                                             if (!(is.list(tmp))) rowSums(tmp) else 0
            sem_object$scale[k] <- .rsf/sqrt(mean(tmp2^2))
        }
        X[,k] <- sem_object$scale[k] * (sem_object$mu[k] + 
            if (!(is.list(tmp))) rowSums(tmp) else 0)  + X[,k] * sem_object$e_var[k]
    }
    if (!is.null(seed_))
        .GlobalEnv$.Random.seed <- seed.bak
    return( if (.rescale) sem_object else X)
}