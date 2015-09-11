
context("simple example")

#library(CAM);library(Matrix) # in case you want to test by hand
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
test_that("basic", { set.seed(1)
    expect_equal(trueDAG, CAM(X, scoreName = "SEMGAM", variableSel = TRUE, pruning = TRUE)$Adj)
})

test_that("order fixation", { set.seed(1)
    expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "force_edge", pruning=T)$Adj)
    expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj)
    expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "force_edge", pruning=T)$Adj))
    expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "emulate_edge", pruning=T)$Adj))
})

test_that("fit model given DAG", { set.seed(1)
    estDAG <- CAM(X, variableSel = TRUE, pruning = TRUE)
    cam1 <- cam.fit(X, estDAG$Adj) 
    expect_less_than(logLikScore(predict(cam1,X)), estDAG$Score+0.00000001)
    
    estDAG <- CAM(X, variableSel = FALSE, pruning = FALSE)
    cam2 <- cam.fit(X, estDAG$Adj) 
    expect_equal(estDAG$Score, logLikScore(predict(cam2,X)))
})

test_that("var test positive", { set.seed(1)
    estDAG_F <- CAM(X,fixedOrders = c(2,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    estDAG_R <- CAM(X,fixedOrders = c(1,2),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    
    cam_F <- cam.fit(X, estDAG_F) 
    cam_R <- cam.fit(X, estDAG_R)
    
    expect_less_than(var.test(cam_F, cam_R)$p.value, 0.05)
})

test_that("var test negative", { set.seed(1)
    estDAG_F <- CAM(X,fixedOrders = c(3,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    estDAG_R <- CAM(X,fixedOrders = c(1,3),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    
    cam_F <- cam.fit(X, estDAG_F) 
    cam_R <- cam.fit(X, estDAG_R)
    
    expect_more_than(var.test(cam_F, cam_R)$p.value, 0.05)
})

test_that("slimming works", { set.seed(1)
    cam <- cam.fit(X, trueDAG) 
    expect_equal(predict(cam, X)$fitted.values, predict(slim(cam), X)$fitted.values)
})

test_that("predicting works", {
    cam <- cam.fit(X, trueDAG)
    cam_fits <- cam$fitted.values
    predicted_fits <- predict(cam, X)$fitted.values
    expect_equal(cam_fits, predicted_fits)
})

test_that("dag to causal order to dag works", {
    causalOrder <- dagToCausalOrder(trueDAG)
    adjacency <- causalOrderToAdjacency(causalOrder)
    expect_true(all(adjacency[as.matrix(trueDAG)]))
})


test_that("", {
    pathMatrixShould <- structure(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, 
                                    TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE), .Dim = c(4L, 4L))
    expect_equal(getPathMatrix(as.matrix(trueDAG)), pathMatrixShould)
    expect_equal(matrix(F, nrow=4, ncol=4), matrix(F, nrow=4, ncol=4))
    expect_equal(matrix(F, nrow=0, ncol=0), matrix(F, nrow=0, ncol=0))
})

