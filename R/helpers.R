dagToCausalOrder <- function(causalDAG)
{
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

causalOrderToAdjacency <- function(causalOrder)
{
    p <- length(causalOrder)
    adjacency <- matrix(FALSE,p,p)
    adjacency[row(adjacency) < col(adjacency)] <- TRUE
    adjacency[causalOrder,causalOrder]
}