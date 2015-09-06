test_that("simple example dataset", {
    #library(CAM)
    #library(Matrix)
  set.seed(1)
  n <- 500
  eps1<-rnorm(n)
  eps2<-rnorm(n)
  eps3<-rnorm(n)
  eps4<-rnorm(n)
  
  x2 <- 0.5*eps2
  x1 <- 0.9*sign(x2)*(abs(x2)^(0.5))+0.5*eps1
  x3 <- 0.8*x2^2+0.5*eps3
  x4 <- -0.9*sin(x3) - abs(x1) + 0.5*eps4
  
  X <- cbind(x1,x2,x3,x4)
  
  trueDAG <- sparseMatrix(i=c(3, 2, 2, 1),
                          j=c(4, 3, 1, 4),dims=c(4,4)) 
  ## x4 <- x3 <- x2 -> x1 
  ##  ^               /
  ##   \_____________/
  ## adjacency matrix:
  ## 0 0 0 1
  ## 1 0 1 0
  ## 0 0 0 1
  ## 0 0 0 0
  
  estDAG <- CAM(X, scoreName = "SEMGAM", numCores = 1, output = F, variableSel = T, 
                pruning = TRUE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001))
  
  expect_equal(trueDAG,estDAG$Adj)
  
  expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "force_edge", pruning=T)$Adj)
  
  expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj)
  
  expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "force_edge", pruning=T)$Adj))
  
  expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "emulate_edge", pruning=T)$Adj))
})
