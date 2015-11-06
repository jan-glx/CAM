#' @export
fitAllOrders <- function(X, nodeModelName = NULL, nodeModelPars = NULL, verbose = FALSE) {
    p <- ncol(X)
    orders <- allOrderPermutations(p)
    results <- data.table(causalOrder = orders,
                          cam = lapply(orders, function(causalOrder) {
                              cam.fit(X, causalDAG=causalOrderToAdjacency(causalOrder), 
                                      nodeModelName = nodeModelName, 
                                      nodeModelPars = nodeModelPars,
                                      verbose = verbose)
                          }))
    results[, ':='(score = sapply(cam, logLik),
                   lli = lapply(cam, logLik_single.cam)
                   )]
    results[, id:=.I]
    setorderv(results, "score", order=-1L)
    
    constraints <- combn(p, 2)
    constraints <- cbind(constraints,constraints[c(2,1),])
    causalOrders = simplify2array(results[, causalOrder])
    mll <- data.table(ii = constraints[1,],
                      jj = constraints[2,],
                      id = apply(constraints, 2, function(constrain) {
                          results[ causalOrders[constrain[2],]<causalOrders[constrain[1],], id[1]]# first one is the one with max logLik
                      }))
    mll <-results[mll,.(causalOrder=causalOrder, score=score, lli=lli, id=id,  ii=ii, jj=jj), on="id"]
    orderIDs2check <- unique(mll[,id])
    if (verbose) {
        cat("of the", ncol(causalOrders), "posible causal orders", length(orderIDs2check), 
            "are optimal under at least one of the", ncol(constraints), "constraints")
        cat("."); flush.console()
    }	
    return(list(fits = results, bestFits = results[id %in% orderIDs2check], mll=mll))
}

mlloc <- function(results, p) {
    constraints <- combn(p, 2)
    constraints <- cbind(constraints,constraints[c(2,1),])
    causalOrders = simplify2array(results[, causalOrder])
    mll <- data.table(ii = constraints[1,],
                      jj = constraints[2,],
                      id = apply(constraints, 2, function(constrain) {
                          results[causalOrders[constrain[2],]<causalOrders[constrain[1],], 
                               id[1]]# first one is the one with max logLik
                      }))
    mll[results, ':='(score=score, lli=lli), on="id"]
    return(mll)
}

