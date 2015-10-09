quick <- Sys.getenv("USERNAME") != "jan" && FALSE

are_causal_orders_compatible_to <- function(trueDAG) {
    s_trueDAG <- substitute(trueDAG)
    function(estDAG) {
        s_estDAG <- substitute(estDAG)
        expectation(areAllCausalOrdersCompatible(estDAG, trueDAG), 
                    paste("there is a causal order of ", s_estDAG, "\n",
                          paste(capture.output(estDAG), collapse="\n"),
                          "\n not compatible with ", s_trueDAG, "\n",
                          paste(capture.output(trueDAG), collapse="\n")
                    ),
                    paste("all orders of", s_estDAG, " are compatible with",  s_trueDAG)
        )
    }
}

expect_causalOrders_compatible_to <- function(object, trueDAG, info = NULL, label = NULL) {
    if (is.null(label)) {
        label <- testthat:::find_expr("object")
    }
    expect_that(object, are_causal_orders_compatible_to(trueDAG))
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
        expect_equal(trueDAG, CAM(X, nodeModelName = nodeModelName, variableSel = TRUE, pruning = TRUE)$Adj)
    })
}
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


# simpler DAG for boostrap tests -------------------------------------------------------------------
context("bootstrap")
p=3
trueDAG <- matrix(FALSE,ncol=p, nrow=p)
trueDAG[matrix(c(1,2,
                 3,3),ncol=2)] <- TRUE
obj <- random_additive_polynomial_SEM(trueDAG, seed_=2)
X <- CAM::simulate_additive_SEM(obj, n=400,seed_ = 3)

test_that("bootstrap test two sided V(3)", {
    if (quick) skip("quick")
    boot_res <- bootstrap.cam(X, matrix(c(2,1), ncol=2),B=100, method = "two-sided") 
    expect_more_than(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam(X, matrix(c(3,1), ncol=2),B=100, method = "two-sided")
    expect_less_than(boot_res$pvalue, 0.05)
})

test_that("bootstrap test one-sided V(3)", {
    if (quick) skip("quick")
    skip()
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(3,1),ncol=2))
    expect_more_than(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,2),ncol=2))
    expect_more_than(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,3),ncol=2))
    expect_less_than(boot_res$pvalue, 0.05)
})

test_that("bootstrap test one-sided V(3) lvl0", {
    if (quick) skip("quick")
    skip()
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(3,1),ncol=2), bs_lvl0 = TRUE)
    expect_more_than(boot_res$pvalue, 0.05)
    boot_res <- bootstrap.cam.one_sided(X, ij = matrix(c(1,3),ncol=2), bs_lvl0 = TRUE)
    expect_less_than(boot_res$pvalue, 0.05)
})
