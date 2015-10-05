#' Calculate a causal order of a DAG \code{causalDAG}
#' 
#' @param causalDAG Adjacency matrix of the DAG whose order is to be computed.
#' @return \eqn{\pi}(\code{causalDAG})
#' @export 
dagToCausalOrder <- function(causalDAG) {
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
    order(causalOrder)
}

#' Calculate adjacency matrix of the fully connected DAG corresponiding to a causal order
#' @param causalOrder Vector of length \eqn{p} where the \eqn{i}'th element is \eqn{\pi(i)}.
#' @export 
causalOrderToAdjacency <- function(causalOrder) {
    p <- length(causalOrder)
    return(outer(causalOrder, causalOrder, '<'))
}


#' Check if there is causal ordering compatible with both of two DAGs
#' \eqn{\exists \pi \in \Pi(G) \cup \Pi(H)}
#' @param G, H Adjacenzy matrices of the two DAGs.
#' @export
#' @return TRUE if there is causal a ordering compatible with both DAGs
existsCompatibleCausalOrder <- function(G, H) {
  !any(xor(getPathMatrix(G),t(getPathMatrix(H))))
}

#' Check if all causal orderings of a DAG \code{G} compatible with DAG \code{H}
#' \eqn{\pi \in \Pi(H) \forall \pi \in \Pi(G)}
#' @param G, H Adjacenzy matrices of the two DAGs.
#' @export
#' @return TRUE if there is causal a ordering compatible with both DAGs
areAllCausalOrdersCompatible <- function(G, H) {
    all(getPathMatrix(G)[as.matrix(H)])
}

#' Get compatible causal order
#' @param estDAG Adjacenzy matrix of the DAG whichs order is to be checked.
#' @param trueDAG Adjacenzy matrix of the DAG which defines the set of "correct" orderings.
#' @export
#' @return TRUE if the causal ordering matches the constrains, otherwise FALSE
getCompatibleCausalOrder <- function(estDAG,trueDAG) {
    dagToCausalOrder(as.matrix(G)|t(as.matrix(H)))
}

#' compute adjacency matrix form edges
#' @export
edges2adj <- function(i, j, p=max(c(i,j))) {
    adj <- matrix(FALSE, p, p)
    adj[cbind(i, j)] <- TRUE
    return(adj)
}

#' Turn a adjacency matrix into a matrix of edges
#' @note equivalent to which(adj, arr.ind = T)
#' @param adj adjacency matrix
#' @return m x 2 matix where each row r denotes an edge from r[,1] to r[,2] and m is the number of edges in adj
#' @export
adj2edges <- function(adj) which(adj, arr.ind = T)


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

#' Compute all causal orderings of a number of nodes
#' @param p Number of nodes that are to be permutated.
#' @return List of order vectors (\eqn{pi})
#' @importFrom gtools permutations
#' @export
allOrderPermutations <- function(p) {
    parray <- permutations(p, p)
    return(split(parray,seq_len(dim(parray)[1])))
}
    
    
    
    
    