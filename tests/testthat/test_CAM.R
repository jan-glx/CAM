quick <- Sys.getenv("USERNAME") != "jan" && FALSE

library(testthat)

expect_causalOrders_compatible_to <- function(object, trueDAG, info = NULL, label = NULL) {
    act <- quasi_label(rlang::enquo(object), label = label, arg = "object")
    tru <- quasi_label(rlang::enquo(trueDAG), arg = "trueDAG")

    expect(areAllCausalOrdersCompatible(act$val, tru$val), 
        paste(
            "there is a causal order of ",  act$lab, "\n",
            paste(capture.output(act$val), collapse="\n"),
            "\n not compatible with ", tru$lab, "\n",
            paste(capture.output(tru$val), collapse="\n")
        )
    )
    invisible(act$val)
}

context("simple example")

#   library(CAM) # in case you want to test by hand
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

trueDAG <- edges2adj(i=c(3, 2, 2, 1),
                     j=c(4, 3, 1, 4)) 
## x4 <- x3 <- x2 -> x1 
##  ^               /
##   \_____________/
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0
for (nodeModelName in c("gam", "poly", "lmboost")){
    test_that(paste0("basic:",nodeModelName), { set.seed(1)
        expect_causalOrders_compatible_to(CAM(X, nodeModelName = nodeModelName)$Adj, trueDAG)
    })
}
test_that("linear runs without errors:", { set.seed(1)
    CAM(X, nodeModelName = "linear")
})

for (nodeModelName in c("gam", "poly", "lmboost")){
    test_that(paste0("with Pruning and variable selection:",nodeModelName), { set.seed(1)
        expect_output(CAM(X, nodeModelName = nodeModelName, variableSel = TRUE, pruning = TRUE, verbose=TRUE), regexp=".+")
        expect_equal(trueDAG, CAM(X, nodeModelName = nodeModelName, variableSel = TRUE, pruning = TRUE)$Adj)
    })
}

test_that("lasso PNS runs without errors:", { set.seed(1)
    CAM(X, pnsMethod = "lasso", variableSel = TRUE)
})

test_that("linear runs without errors with pruning and PNS:", { set.seed(1)
    CAM(X, nodeModelName = "linear", variableSel = TRUE, pruning = TRUE)
})

test_that("order fixation", { set.seed(1)
    expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "force_edge", pruning=TRUE)$Adj)
    expect_equal(trueDAG,CAM(X,fixedOrders = c(2,1),orderFixationMethod = "emulate_edge", pruning=TRUE)$Adj)
    expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "force_edge", pruning=TRUE)$Adj))
    expect_true(!all(trueDAG==CAM(X,fixedOrders = c(1,2),orderFixationMethod = "emulate_edge", pruning=TRUE)$Adj))
})

test_that("logLik CAM equal to logLik cam.fit", { set.seed(1)
    estDAG <- CAM(X, variableSel = TRUE, pruning = FALSE)
    cam1 <- cam.fit(X, estDAG$Adj) 
    expect_equal(as.numeric(logLik(cam1)), estDAG$score)
    expect_equal(as.numeric(logLik(predict(cam1, X))), estDAG$score)
})

test_that("cam.fit() equal to predict(cam.fit())", { set.seed(1)
    cam1 <- cam.fit(X, trueDAG) 
    expect_equal(as.numeric(logLik(predict(cam1, X))), as.numeric(logLik(cam1)))
})

test_that("var test positive", { set.seed(1)
    estDAG_F <- CAM(X,fixedOrders = c(2,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    estDAG_R <- CAM(X,fixedOrders = c(1,2),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    
    cam_F <- cam.fit(X, estDAG_F) 
    cam_R <- cam.fit(X, estDAG_R)
    
    expect_lt(var.test(cam_F, cam_R)$p.value, 0.05)
})

test_that("var test negative", { set.seed(1)
    estDAG_F <- CAM(X,fixedOrders = c(3,1),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    estDAG_R <- CAM(X,fixedOrders = c(1,3),orderFixationMethod = "emulate_edge", pruning=T)$Adj
    
    cam_F <- cam.fit(X, estDAG_F) 
    cam_R <- cam.fit(X, estDAG_R)
    
    expect_gt(var.test(cam_F, cam_R)$p.value, 0.05)
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


# simpler DAG for boostrap tests -------------------------------------------------------------------
context("bootstrap")
p=3
trueDAG <- matrix(FALSE,ncol=p, nrow=p)
trueDAG[matrix(c(1,2,
                 3,3),ncol=2)] <- TRUE
sem_object <- random_additive_polynomial_SEM(trueDAG, seed_=5)
sem_object <- rescale_sem_object(sem_object, seed_ = 4)
X <- CAM::simulate_additive_SEM(sem_object, n=400,seed_ = 3)
pairs(X)

test_that("bootstrap test two sided V(3)", {
    if (quick) skip("quick")
    boot_res <- bootstrap.cam(X, matrix(c(2,1), ncol=2),B=100, method = "two-sided") 
    expect_gt(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam(X, matrix(c(3,1), ncol=2),B=100, method = "two-sided")
    expect_lt(boot_res$pvalue, 0.05)
})

test_that("bootstrap test one-sided V(3)", {
    if (quick) skip("quick")
    skip("not implemented anymore")
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(3,1),ncol=2))
    expect_gt(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,2),ncol=2))
    expect_gt(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,3),ncol=2))
    expect_lt(boot_res$pvalue, 0.05)
})

test_that("bootstrap test one-sided V(3) lvl0", {
    if (quick) skip("quick")
    skip("not implemented anymore")
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(3,1),ncol=2), bs_lvl0 = TRUE)
    expect_gt(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,3),ncol=2), bs_lvl0 = TRUE)
    expect_lt(boot_res$pvalue, 0.05)
})

test_that("bootstrap test one-sided V(3) lvl0", {
    if (quick) skip("quick")
    skip("not implemented anymore")
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(3,1),ncol=2), bs_lvl0 = TRUE)
    expect_gt(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,3),ncol=2), bs_lvl0 = TRUE)
    expect_lt(boot_res$pvalue, 0.05)
})

test_that("direct bootstrap test one-sided V(3) works", {
    if (quick) skip("quick")
    dt <- bootstrap.cam.one_sided_direct(X)[[1]]
    dt[, true:=TRUE]
    dt[(ii==1 & jj==3)|(ii==2 & jj==3), true:=FALSE]
    expect_true(all(dt[, xor(p.value<=0.05,true)]))
})

context("edge cases")
set.seed(1)
n <- 19
eps1<-rnorm(n)
eps2<-rnorm(n)
eps3<-rnorm(n)
eps4<-rnorm(n)

x2 <- 0.5*eps2
x1 <- 0.9*sign(x2)*(abs(x2)^(0.5))+0.5*eps1
x3 <- 0.8*x2^2+0.5*eps3
x4 <- -0.9*sin(x3) - abs(x1) + 0.5*eps4

X <- cbind(x1,x2,x3,x4)

trueDAG <- edges2adj(i=c(3, 2, 2, 1),
                     j=c(4, 3, 1, 4)) 
## x4 <- x3 <- x2 -> x1 
##  ^               /
##   \_____________/
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0
for (nodeModelName in c("gam", "poly", "lmboost")){
    test_that(paste0("low n - basic:",nodeModelName), { set.seed(1)
        expect_true({CAM(X, nodeModelName = nodeModelName);TRUE})
    })
}