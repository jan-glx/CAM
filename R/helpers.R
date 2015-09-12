#' @export 
dagToCausalOrder <- function(causalDAG)
{
    causalDAG <- as.matrix(causalDAG)
    p <- nrow(causalDAG)
    idx <- 1:p
    causalOrder <- integer()
    while (length(causalOrder) <p) {
        sink_idx <- which(rowSums(causalDAG)==0)
        causalOrder <- c(idx[sink_idx], causalOrder)  #slow but should be fine for small/sparse graphs
        causalDAG <- causalDAG[-sink_idx,-sink_idx, drop=F]
        idx <- idx[-sink_idx]
    }
    causalOrder
}

#' @export 
causalOrderToAdjacency <- function(causalOrder)
{
    p <- length(causalOrder)
    adjacency <- matrix(FALSE,p,p)
    adjacency[row(adjacency) < col(adjacency)] <- TRUE
    adjacency[causalOrder,causalOrder]
}


#' Check if the causal ordering of an DAG (estDAG) matches the ordering of an other DAG (trueDAG) 
#' using their adjacency matrices
#' @param estDAG Adjacenzy matrix of the DAG whichs order is to be checked
#' @param trueDAG Adjacenzy matrix of the DAG which defines the set of "correct" orderings
#' @export
#' @return TRUE if the causal ordering matches the constrains, otherwise FALSE
checkCausalOrder <- function(estDAG,trueDAG)
{
  all(as.matrix(estDAG)[as.matrix(trueDAG)])
}



#' Compute the path matrix of a graph G from the adjacency matrix of G
#' @param adj adjacency matrix to compute the pathmatrix for. Must be a \eqn{p\times{}p} matrix that
#'   is TRUE at the edges of G and FALSE else where
#' @export
getPathMatrix <- function(adj)
{
  path_matrix=matrix(FALSE, nrow(adj), ncol(adj) )
  diag(path_matrix) <- TRUE
  edges=which(adj,arr.ind = T)
  for(i in seq_along(edges[,1])){
    DE_of_and_child <- which(path_matrix[edges[i,2],])
    AC_of_and_parent <- which(path_matrix[,edges[i,1]])
    path_matrix[AC_of_and_parent,DE_of_and_child] <- TRUE
  }
  return(path_matrix)
}