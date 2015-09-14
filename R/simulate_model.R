#' @export
random_additive_polynomial_SEM <- function(trueDAG, degree=3, seed_=1){
  p <- ncol(trueDAG)

  set.seed(seed_)
  rand_poly <- function(degree){
    coef <- rnorm(degree+1, 0, 1)
    f<- function(x){ t(sapply(x, `^`, (0:degree))) %*% coef}
    attr(f, "coef") <- coef
    f
  }

  f_jk <- matrix(list(),p,p)
  f_jk[trueDAG>0] <- lapply(rep(degree,sum(trueDAG)), FUN=rand_poly)
  list(trueDAG=trueDAG, f_jk=f_jk, p=p)
}

#' @export
simulate_additive_SEM <- function(sem_object, n=500, scaling=FALSE, seed_=1){
  set.seed(seed_)

  p <- sem_object$p
  trueDAG <- sem_object$trueDAG
  f_jk <- sem_object$f_jk

  hsplit <- function(x){
    x = as.matrix(x)
    split(x, col(x))
  }

  X <- matrix(rnorm(n*p), ncol=p) # noise
  for(k in dagToCausalOrder(as.matrix(sem_object$trueDAG))) {
    tmp <- mapply(do.call,
               f_jk[trueDAG[,k]>0,k],
               lapply(hsplit(X[,trueDAG[,k]>0]),list))
    if (!(is.list(tmp))) {
      X[,k] <- X[,k] + rowSums(if(scaling) scale(tmp) else tmp)
    }
  }
  X
}