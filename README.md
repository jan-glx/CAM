# Causal Additive Model (CAM) [![Build Status](https://travis-ci.org/jan-glx/CAM.svg)](https://travis-ci.org/jan-glx/CAM) [![codecov.io](http://codecov.io/github/jan-glx/CAM/coverage.svg?branch=master)](http://codecov.io/github/jan-glx/CAM?branch=master)
**And R package to infer the causal DAG under the assumtion of an additive model.**

The code takes an n x p data matrix and fits a Causal Additive Model (CAM) for estimating the causal structure of the underlying process. The output is a p x p adjacency matrix (a one in entry (i,j) indicates an edge from i to j). 

Details of the algorithm can be found in: P. Bühlmann, J. Peters, J. Ernest: "CAM: Causal Additive Models, high-dimensional order search and penalized regression", Annals of Statistics 42:2526-2556, 2014.
## Changes
* adds bootstrap based hypothesis test for causal effects
* adds functions to simulate related models
* adds some documentation & unit-tests
* simplifies the code base & speeds up model fits (using data.table)

## Installation
```r
# install.packages("remotes")
remotes::install_github("jan-glx/CAM")
```

## Usage
```
library(CAM)
```
### simple model fit
#### # set-up data
```
set.seed(1)
n <- 5000
eps1 <- rnorm(n)
eps2 <- rnorm(n)
eps3 <- rnorm(n)
eps4 <- rnorm(n)

x2 <- 0.5 * eps2
x1 <- 0.9 * sign(x2) * (abs(x2) ^ (0.5)) + 0.5 * eps1
x3 <- 0.8 * x2 ^ 2 + 0.5 * eps3
x4 <- -0.9 * sin(x3) - abs(x1) + 0.5 * eps4

X <- cbind(x1, x2, x3, x4)

trueDAG <- edges2adj(i = c(3, 2, 2, 1),
                     j = c(4, 3, 1, 4)) 
## x4 <- x3 <- x2 -> x1 
##  ^               /
##   \_____________/
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0
```
#### # fit CAM
```
fit1 <- CAM(X)
fit1$Adj # fitted DAG
#>       [,1]  [,2]  [,3]  [,4]
#> [1,] FALSE FALSE FALSE  TRUE
#> [2,]  TRUE FALSE  TRUE  TRUE
#> [3,]  TRUE FALSE FALSE  TRUE
#> [4,] FALSE FALSE FALSE FALSE

# check if causal ordering of inferred causal DAG is compatible with true causal DAG
areAllCausalOrdersCompatible(fit1$Adj, trueDAG)
#> [1] TRUE
```
### Bootstrap-based hypothesis test
#### # set-up example DAG
```
p=3
trueDAG <- matrix(FALSE, ncol=p, nrow=p)
trueDAG[matrix(c(1, 2,
                 3, 3), ncol=2)] <- TRUE
## 1 -> 3 <- 2
trueDAG
#>       [,1]  [,2]  [,3]
#> [1,] FALSE FALSE  TRUE
#> [2,] FALSE FALSE  TRUE
#> [3,] FALSE FALSE FALSE
```
#### # simulate example data
```
sem_object <- random_additive_polynomial_SEM(trueDAG, seed_ = 5)
sem_object <- rescale_sem_object(sem_object, seed_ = 4)
X <- simulate_additive_SEM(sem_object, n = 400, seed_ = 3)
pairs(X)
```
#### # perform tests
```
# Test H0: 2 <~/~> 1 (TRUE)
bootstrap.cam(X, matrix(c(2, 1), ncol=2), B=100, method = "two-sided")$pvalue 
#> [1] 0.86

# Test H0: 2 <~/~> 3 (FALSE)
bootstrap.cam(X, matrix(c(2, 3), ncol=2), B=100, method = "two-sided")$pvalue 
#> [1] 0
```

## Authors
This basically the package from CRAN by J. Peters, J. Ernest with some minor changes by [@jan-glx](https://github.com/jan-glx) and [@nignatiadis](https://github.com/nignatiadis).
