

#' Causal Additive Model
#' 
#' fits a causal additive model using the CAM algorithm, see references below
#' 
#' The code fits a CAM model. See the references below for more details.
#' Identifiability results for the model class can be found in
#' 
#' J. Peters, J. Mooij, D. Janzing, B. Sch\"olkopf: Causal Discovery with
#' Continuous Additive Noise Models, JMLR 15:2009-2053, 2014.
#' 
#' @param X nxp matrix of training inputs (n data points, p dimensions)
#' @param scoreName specifies the model type which is used to compute the
#' score. Default is "SEMGAM" which assumes a generalized additive model class.
#' Other options include "SEMLIN" which fits a linear model.
#' @param parsScore additional parameters can be supported to the score
#' function.
#' @param numCores specifies the number of cores that can be used for
#' computation.
#' @param maxNumParents specifies the maximal number of parents that are
#' allowed in the model.
#' @param verbose shall verbose be printed to the console (TRUE/FALSE)
#' @param variableSel specifies whether initial variable selection (Step 1 of
#' CAM algorithm) shall be performed (TRUE) or not (FALSE). Initial variable
#' selection reduces the number of possible parents for a given node and
#' therefore enables computing the causal structure for large p.
#' @param variableSelMethod specifies the method that is used for variable
#' selection. Default is selGamBoost which uses the gamboost function from
#' mboost package. Other options include: selGam (gam() from mgcv), selLm based
#' on linear regression, selLasso based on Lasso regression from package
#' glmnet.
#' @param variableSelMethodPars optional parameters to modify settings of the
#' selection method.
#' @param pruning specifies whether pruning (Step 3 of CAM algorithm) shall be
#' performed (TRUE) or not (FALSE). Pruning reduces the number of edges in the
#' estimated causal structure.
#' @param pruneMethod specifies the method used for the pruning step. Default
#' is selGAM which is based on the gam() function from the mgcv package.
#' @param pruneMethodPars optional parameters to tune the pruning step.
#' @param intervData boolean that indicates whether we use interventional data.
#' @param intervMat the matrix intervMat has the same dimension as X. entry
#' (i,j) == TRUE indicates that in experiment i, variable j has been intervened
#' on.
#' @return list of attributes of the final estimated causal structure
#' \item{Adj}{adjacency matrix of estimated causal graph} \item{Score}{Total
#' edge score of estimated graph} \item{timesVec}{Vector containing various
#' time measurements for execution times of the individual steps of the CAM
#' algorithm}
#' @author Jonas Peters <jonas.peters@@tuebingen.mpg.de> and Jan Ernest
#' <ernest@@stat.math.ethz.ch>
#' @references P. B\"uhlmann, J. Peters, J. Ernest: CAM: Causal Additive
#' Models, high-dimensional Order Search and Penalized Regression Annals of
#' Statistics 42:2526-2556, 2014.
#' @keywords Causality Regression Additive Noise Models Restricted Structural
#' Equation Model
#' @examples
#' 
#' n <- 500
#' eps1<-rnorm(n)
#' eps2<-rnorm(n)
#' eps3<-rnorm(n)
#' eps4<-rnorm(n)
#' 
#' x2 <- 0.5*eps2
#' x1 <- 0.9*sign(x2)*(abs(x2)^(0.5))+0.5*eps1
#' x3 <- 0.8*x2^2+0.5*eps3
#' x4 <- -0.9*sin(x3) - abs(x1) + 0.5*eps4
#' 
#' X <- cbind(x1,x2,x3,x4)
#' 
#' trueDAG <- cbind(c(0,1,0,0),c(0,0,0,0),c(0,1,0,0),c(1,0,1,0))
#' ## x4 <- x3 <- x2 -> x1 -> x4
#' ## adjacency matrix:
#' ## 0 0 0 1
#' ## 1 0 1 0
#' ## 0 0 0 1
#' ## 0 0 0 0
#' 
#' estDAG <- CAM(X, pruning = TRUE)
#' 
#' cat("true DAG:\n")
#' show(trueDAG)
#' 
#' cat("estimated DAG:\n")
#' show(estDAG$Adj)
#' 
#' @export CAM
#' @import data.table
CAM <-
    function(X, nodeModelName = c("gam", "lasso", "poly", "linear", "lmboost"), scoreFunction = logLik, 
             nodeModelPars = NULL, maxNumParents = ncol(X) - 1, #min(dim(X)[2] - 1, round(dim(X)[1]/20)),
             fixedOrders = matrix(ncol = 2, nrow = 0),
             orderFixationMethod = c("force_edge", "emulate_edge", "free"), 
             variableSel = FALSE, pnsMethod = NULL, pnsModelPars = NULL, pnsSelectionPars = NULL,
             pruning = FALSE, pruneMethod = NULL, pruneModelPars = NULL, pruneSelectionPars = NULL,
             intervData = FALSE, intervMat = NULL,
             numCores = 1, verbose = FALSE) {

    # input checking -------------------------------------------------------------------------------
    if (!is.function(nodeModelName)) nodeModelName <- match.arg(nodeModelName)
    if (is.null(nodeModelName)) nodeModelName <- "gam"
    orderFixationMethod <- match.arg(orderFixationMethod)
    if (class(fixedOrders) == "numeric" && length(fixedOrders) == 2){
      fixedOrders <- matrix(as.integer(fixedOrders), nrow=1)
    } else if (is.null(fixedOrders)) {
      fixedOrders <- matrix(ncol=2,nrow=0)
    } else if (!is.matrix(fixedOrders) || ncol(fixedOrders) != 2){
      stop("fixedOrders must be a [nrow X 2] matrix where pi(fixedOrders[i,1]) < pi(fixedOrders[i,2]), and pi is a causal ordering of the nodes")
    }
  
    if(verbose) cat("number of cores:", numCores, "\n")
    if (is.list(X)||is.data.frame(X)||is.data.table(X)) setDT(X)
    else X <- as.data.table(X)
    p <- dim(X)[2]
    
    ## STEP 1: variable selection
    if(variableSel) {
        selMat <- preliminaryNeighborSel(X, method = pnsMethod, modelPars = pnsModelPars, 
                                         selectionPars = pnsSelectionPars, intervData = intervData, 
                                         intervMat = intervMat, numCores = numCores, verbose = verbose)
    } else {
        selMat <- matrix(TRUE, p,p)
    }

    # STEP 2: Include Edges ------------------------------------------------------------------------
    cam <- incEdge(X, nodeModelName = nodeModelName, nodeModelPars = nodeModelPars,
                   scoreFunction = scoreFunction, maxNumParents = maxNumParents,
                   fixedOrders = fixedOrders, orderFixationMethod = orderFixationMethod, 
                   selMat=selMat, intervData = intervData, intervMat = intervMat, 
                   numCores = numCores, verbose = verbose)
    
    # STEP 3: Prune the DAG ------------------------------------------------------------------------
    if(pruning) cam$Adj <- prune(X=X,G=cam$Adj, method = pruneMethod, modelPars = pruneModelPars, 
                                 selectionPars = pruneSelectionPars, intervData = intervData, 
                                 intervMatrix = intervMatrix, numCores = numCores, verbose = verbose)  
    
    if(verbose) cat("\n ... done. \n ")
    return(c(cam, list(score = sum(cam$scoreNodes))))  
}

prune <- function(X, G, method = NULL, modelPars = NULL, selectionPars = NULL, intervData = FALSE, 
                  intervMatrix = NULL, numCores = 1, verbose = FALSE) {
    p <- ncol(X)
    if (is.null(method)) method <- "gam"
    if (verbose) cat("\n Performing pruning ... \n ")
    if (verbose && intervData) cat("The pruning is done with the observational data only.\n")
    X2 <- if(intervData) X[rowSums(intervMat) == 0,] else X
    
    finalG <- matrix(FALSE, p, p)
    for(j in 1:p) {
        parents <- which(G[,j])
        if (length(parents) > 0) {
            if (verbose) cat("pruning variable:", j, "\n", "considered parents:", parents, "\n")
            selectedPar <- selectParents(X2, j, parents_of_j = parents, method = method, 
                                         modelPars = modelPars, selectionPars = selectionPars, verbose = verbose)
            finalG[selectedPar, j] <- TRUE
        }
    }
    return(finalG)
}

#' A matrix selMat is constructed. Entry (i,j) being one means that i is a possible parent of j.
preliminaryNeighborSel <- function(X, method = NULL, modelPars = NULL, selectionPars = NULL, 
                                   intervData = FALSE, intervMat = NULL, numCores = 1, verbose = F) {
    p <- ncol(X)
    if (is.null(method)) method <- "mboost"
    if (verbose && intervData) cat("The preliminary neighbourhood selection is done with the observational data only.\n")
    X2 <- if(intervData) X[rowSums(intervMat) == 0,] else X # we should drop rows for each j individually here
    
    moreArgs <- list(X = X2, method = method, modelPars = modelPars, selectionPars = selectionPars, verbose = verbose)
    if(numCores == 1) {
        selMat <-             mapply(selectParents, MoreArgs = moreArgs, j=1:p)
    } else {
        selMat <- parallel::mcmapply(selectParents, MoreArgs = moreArgs, j=1:p, mc.cores = numCores)
    }
    # The next line includes j as a possible parent of i if i is considered a possible parent of j
    # selMat <- selMat | t(selMat)
    if(verbose) {
        cat("Instead of p2^(p-1) -Sillander- ",p*2^(p-1) ," we have ", sum(2^rowSums(selMat)), "\n")
        cat("Greedy, on the other hand, is computing ",sum(selMat) ," entries. \n")
        if(p<30) {
            cat("This is the matrix of possible parents after the first step.\n")
            show(selMat)
        }
        cat("Object size of selmat: ", object.size(selMat), "\n")
    }
    return(selMat)
}

selectParents <- function(X, j, parents_of_j = setdiff(1:ncol(X),j), method, modelPars = NULL, 
                          selectionPars = NULL, verbose = FALSE) {
    p <- ncol(X)
    if(verbose) cat("Performing variable selection for variable ", j, ": \n")
    selVec <- rep(FALSE, p)
    if(p >= 2) {
        default_pars <- switch(method,
                               mboost  = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),
                               lmboost = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),
                               gam = list(cutOffPVal = 0.001),
                               lm = list(cutOffPVal = 0.001),
                               lasso = list(),
                               list())
        if (is.null(selectionPars)) selectionPars <- list()
        selectionPars <-  modifyList(default_pars, selectionPars)
        nodeModel <- fitNode(X, j, parents_of_j, method = method, pars = modelPars)
        tmp <- do.call(selectFeatures, c(list(object=nodeModel), selectionPars))
        selVec[parents_of_j] <- do.call(selectFeatures, c(list(object=nodeModel), selectionPars))
    }
    return(selVec)
}

selectFeatures <- function(object, ...) UseMethod("selectFeatures")

selectFeatures.gam <- function(object, cutOffPVal, verbose = FALSE) {
    pValVec <- summary(object)$s.pv
    if(verbose) cat("vector of p-values:", pValVec, "\n")
    return(pValVec < cutOffPVal)
}

selectFeatures.lm <- function(object, cutOffPVal, verbose = FALSE) {
    pValVec <- coef(summary(object))[-1,"Pr(>|t|)"] 
    if(verbose) cat("vector of p-values:", pValVec, "\n")
    return(pValVec < cutOffPVal)
}

selectFeatures.glmnet <- function(object, verbose = FALSE) {
    return(object$beta != 0)
}

selectFeatures.mboost <- function(object, atLeastThatMuchSelected, atMostThatManyNeighbors, verbose = FALSE) {
    dt <- data.table(variable = object$xselect())
    dt <- dt[,.N,by=variable][,.(variable, N=N/nrow(dt))][N>0,]
    if(verbose) {
        cat("The following variables \n")
        show(dt$variable)
        cat("... have been selected in the following fractions: \n")
        show(dt$N)
    }
    dt <- dt[N>=atLeastThatMuchSelected]
    dt <- dt[order(N)[0:min(atMostThatManyNeighbors,nrow(dt))]]
    if(verbose) {
        cat("We finally choose as possible parents: \n")
        show(dt$variable)
        cat("\n")
    }
    tmp <- rep(FALSE, length(object$baselearner))
    tmp[dt[,variable]] <- TRUE
    return(tmp)
}

#' @export 
#' @import data.table
fitNode <- function(X, j, parents_of_j, method = c("gam", "lasso", "poly", "linear", "lmboost", "mboost"), 
                    pars = list(), verbose = FALSE) {
    if (is.function(method)) {
        return(do.call(method, c(list(X=X, j=j, parents_of_j=parents_of_j), pars)))
    }
    method <- match.arg(method)
    default_pars <- switch(method,
                           gam = list(numBasisFcts = 10),
                           poly = list(degree=3),
                           list())
    if (is.null(pars))
        pars <- list()
    pars <-  modifyList(default_pars, pars)
    if(method=="linear") {
        method <- "poly"
        pars$degree <- 1
    }
    
    make_additive_formula <- function (summand_expr = "%s") {
        f <- paste0(colnames(X)[j],"~ 1", paste0(sprintf(paste0("+",summand_expr), 
                                                         colnames(X)[parents_of_j]), 
                                                 collapse=""))
        formula(f)
    }
    res <- switch(method,
                  gam = {
                      nobs <- nrow(X)
                      max_numBasisFcts <-ceiling(nobs/(3*length(parents_of_j)))
                      if (max_numBasisFcts < pars$numBasisFcts) {
                          warning("changed number of basis functions from", pars$numBasisFcts,"to    ", 
                              max_numBasisFcts, "    in order to have enough samples per basis function\n")
                          pars$numBasisFcts <- max_numBasisFcts
                      }
                      f <- make_additive_formula(paste0("s(%s, k=", pars$numBasisFcts, ")"))
                      res <- try(mgcv::gam(formula=f, data=X),silent = TRUE)
                      if(typeof(res) == "logical" || inherits(res, "try-error")) {
                          warning("There was some error with gam. The smoothing parameter is set to zero.\n")
                          f <-make_additive_formula(paste0("s(%s, k=", pars$numBasisFcts, ", sp=0)"))
                          res <- mgcv::gam(formula=f, data=X)
                      }
                      res
                  },
                  lasso = {
                      if (length(parents_of_j)>2) {
                          cvres <- glmnet::cv.glmnet(as.matrix(X[, parents_of_j, with=F]), as.matrix(X[, j, with=F]))
                          glmnet::glmnet(as.matrix(X[, parents_of_j, with=F]), as.matrix(X[, j, with=F]), lambda = cvres$lambda.1se)
                      } else {
                          glm(formula = make_additive_formula(), data = X)
                      }
                  },
                  poly = {
                      f <- make_additive_formula(paste0("poly(%s,degree=", pars$degree, ",raw=TRUE)"))
                      lm(formula=f, data=X)
                  },
                  lmboost = {
                      op <- options(warn=-1)
                      on.exit(options(op))
                      mboost::glmboost(make_additive_formula("%s"), X, center = TRUE)
                  },
                  mboost = {
                      bl <- lapply(X[, parents_of_j, with=F], mboost::bbs)
                      mboost::mboost_fit(bl, X[, j, with=F][[1]])
                  },
                  linear = {
                      lm(make_additive_formula("%s"), X)
                  }
    )
    return(res)
}


#' @export  
#' @import data.table
cam.fit <- function(X, causalDAG=NULL, nodeModelName = c("gam", "lasso", "poly", "linear", "lmboost"), 
                    nodeModelPars = NULL, verbose = FALSE)
{
    if (is.null(causalDAG)) stop("Not implemented here. Use CAM(...) instead.")
    if (!is.function(nodeModelName)) nodeModelName <- match.arg(nodeModelName) 
    if (is.list(X)||is.data.frame(X)||is.data.table(X)) setDT(X)
    else X <- as.data.table(X)
    if (!any(class(causalDAG) %in% c("matrix", "ngCMatrix"))) stop("causalDAG must be of class 'matrix' or 'ngCMatrix'")
    p <- nrow(causalDAG)
    nodeModels <- lapply(setNames(as.list(1:p), colnames(X)), function(j) 
        fitNode(X=X, j=j, parents_of_j = which(causalDAG[,j]), method = nodeModelName, 
                pars = nodeModelPars, verbose=verbose))
    vals <- setDT(lapply(nodeModels,fitted.values))
    cam <- list(call= match.call(), causalDAG = causalDAG, p = p, nodeModels = nodeModels, data = X, 
                fitted.values = vals, 
                df.residual = sapply(nodeModels, "[[", "df.residual"),
                est.df = sapply(nodeModels,  function(x) sum(x$est.residual))
                )
    class(cam) <- "cam"
    return(cam)
}



#' @export 
#' @import data.table 
#' @importFrom stats predict
predict.cam <- function(object, newdata, ...)
{
    object$data <- as.data.table(newdata)
    object$fitted.values <- setDT(lapply(object$nodeModels, predict, object$data))
    object$df.residual <- nrow(object$data)
    return(object)
}

#' @export  
#' @importFrom stats residuals
residuals.cam <- function(object, ...) object$data - object$fitted.values

#' @export 
#' @importFrom stats logLik
logLik.cam <- function(object, ...) {
    res <- residuals(object)
    N <- nrow(res)
    #val <-  sum(-1/2*N*(log(sapply(res^2,mean))+log(2*pi)+1))
    logLiks <- lapply(object$nodeModels, logLik)
    val <- do.call(sum, logLiks)
    attr(val, "df") <- sum(sapply(logLiks,function(obj)attr(obj, "df"))) #there are more, but let's ignore them
    class(val) <- "logLik"
    return(val)
}

#' @export 
logLikScore <- function(object, ...) UseMethod("logLikScore")

#' @export 
logLikScore.default <- function(object, ...) -sum(log(var(residuals(object))))
#' @export 
logLikScore.cam <- function(object, ...) -sum(log(sapply(residuals(object), var)))


#' @export 
print.cam <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCausal DAG:\n")
    print(x$causalDAG)
    cat("Log likelihod score:", logLikScore(x))
    invisible(x)
}

#' @export 
var.residuals <- function(object, ...) UseMethod("var.residuals")

#' @export 
var.residuals.cam <- function(object, ...) {
    apply(residuals(object)^2, 2, sum)/object$df.residual
}

# @export
#setClass("cam")

# @export
#setClass("slimmed.cam")


#' @export 
#' @importFrom stats var.test
var.test.cam <- function (x, y, ratio = 1, alternative = c("two.sided", "less", "greater"), 
                          conf.level = 0.95, ...) 
{
    x2 <- list(df.residual = sum(x$df.residual), residuals = sqrt(sum(var.residuals(x))*sum(x$df.residual)))
    y2 <- list(df.residual = sum(y$df.residual), residuals = sqrt(sum(var.residuals(y))*sum(y$df.residual)))
    class(x2) <- "lm"
    class(y2) <- "lm"
    return(stats::var.test(getLmForVarTest(x), getLmForVarTest(y), ratio, alternative, conf.level, ...))
}

#' @export 
getLmForVarTest <- function(x)
{
    x <- list(df.residual = sum(x$df.residual), residuals = sqrt(sum(var.residuals(x))*sum(x$df.residual)))
    class(x) <- "lm"
    return(x)
}

#' @export 
slim <- function(object) UseMethod("slim")

#' @export 
slim.gam <- function(object, ...)
{
    object$offset <- NULL
    object$model <- NULL
    for (i in seq_along(object$smooth)) {
        #object$smooth[[i]]$UZ <- NULL # prediction does not work
        #object$smooth[[i]]$Xu <- NULL # without
    }
    object$y <- NULL
    object$weights <- NULL
    object$prior.weights <- NULL
    object$linear.predictors <- NULL
    object$fitted.values <- NULL
    object$residuals <- NULL
    class(object) <- c("slimmed.gam",class(object))
    return(object)
}

#' @export 
logLikScore.slimmed.cam <- function(object, ...) return(object$logLikScore)

#' @export 
var.residuals.slimmed.cam <- function(object, ...) return(object$var.residuals)

#' @export 
slim.cam <- function(object, ...)
{
    object$var.residuals <- var.residuals(object)
    #object$logLik <- logLik(object)
    object$logLikScore <- logLikScore(object)
    object$fitted.values <- NULL
    object$data <- NULL
    object$nodeModels <- lapply(object$nodeModels, slim)
    class(object) <- c("slimmed.cam",class(object))
    return(object)
}

