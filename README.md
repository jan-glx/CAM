# Causal Additive Model (CAM) [![Build Status](https://travis-ci.org/jan-glx/CAM.svg)](https://travis-ci.org/jan-glx/CAM) [![codecov.io](http://codecov.io/github/jan-glx/CAM/coverage.svg?branch=master)](http://codecov.io/github/jan-glx/CAM?branch=master)
**And R package to infer the causal DAG under the assumtion of an additive model.**

The code takes an n x p data matrix and fits a Causal Additive Model (CAM) for estimating the causal structure of the underlying process. The output is a p x p adjacency matrix (a one in entry (i,j) indicates an edge from i to j). 

Details of the algorithm can be found in: P. Bühlmann, J. Peters, J. Ernest: "CAM: Causal Additive Models, high-dimensional order search and penalized regression", Annals of Statistics 42:2526-2556, 2014.
## Changes
* adds bootstrap based hypothesis test for causal effects
* adds functions to simulate related models
* adds some documentation & unit-tests
* simplifies the code base & speeds up model fits (using data.table)

## Authors
This basically the package from CRAN by J. Peters, J. Ernest with some minor changes by jan-glx and nignatiadis
