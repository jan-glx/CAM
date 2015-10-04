context("helpers")

test_that("dagToCausalOrder", {
    trueDAG <- edges2adj(i=c(3, 3, 1),
                         j=c(1, 2, 2))
    expect_equal(dagToCausalOrder(trueDAG), c(2, 3, 1))
    trueDAG <- edges2adj(i=integer(),
                         j=integer(), p = 1) 
    expect_equal(dagToCausalOrder(trueDAG), 1L)
    
})

test_that("causalOrderToAdjacency", {
    causalOrder <- 1
    expect_equal(causalOrderToAdjacency(causalOrder), matrix(FALSE, 1, 1))
    causalOrder <- c(2, 1)
    expect_equal(causalOrderToAdjacency(causalOrder), edges2adj(i=c(2),
                                                                j=c(1)))
    causalOrder <- c(2, 3, 1)
    expect_equal(causalOrderToAdjacency(causalOrder), edges2adj(i=c(3, 3, 1),
                                                                j=c(1, 2, 2)))
    causalOrder <- c(2, 1)
    expect_equal(causalOrderToAdjacency(causalOrder), edges2adj(i=c(2),
                                                                j=c(1)))
    
})
          

test_that("dagToCausalOrder o causalOrderToAdjacency <=> identity", {
    causalOrder <- sample(10)
    expect_equal(dagToCausalOrder(causalOrderToAdjacency(causalOrder)), causalOrder)
})

trueDAG <- edges2adj(i=c(3, 2, 2, 1),
                     j=c(4, 3, 1, 4))

test_that("dag to causal order to dag works", {
    causalOrder <- dagToCausalOrder(trueDAG)
    fullDAG <- causalOrderToAdjacency(causalOrder)
    expect_true(areAllCausalOrdersCompatible(fullDAG, trueDAG))
})

test_that("path Matrix", {
    pathMatrixShould <- structure(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, 
                                    TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE), .Dim = c(4L, 4L))
    expect_equal(getPathMatrix(as.matrix(trueDAG)), pathMatrixShould)
    expect_equal(matrix(F, nrow=4, ncol=4), matrix(F, nrow=4, ncol=4))
    expect_equal(matrix(F, nrow=0, ncol=0), matrix(F, nrow=0, ncol=0))
})