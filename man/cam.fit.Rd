% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAM.R
\name{cam.fit}
\alias{cam.fit}
\title{Fit a Causal Additive Model (CAM) based on a given DAG
Fit a Causal Additive Model (CAM) based on a given DAG \code{causalDAG} using a specified model 
class (\code{nodeModelName}, \code{nodeModelPars}) for the nodes.}
\usage{
cam.fit(
  X,
  causalDAG = NULL,
  nodeModelName = c("gam", "lasso", "poly", "linear", "lmboost"),
  nodeModelPars = NULL,
  numCores = 1,
  verbose = FALSE
)
}
\arguments{
\item{X}{Data to which the CAM gets fitted. Schould be a data.table but may be a list of numerics
or data.frame.}

\item{causalDAG}{Adjaziency matrix specifying the causal DAG of the CAM to be fitted.}

\item{nodeModelName}{String specifying the model class/ fitting method.}

\item{nodeModelPars}{Additional paramters for the fitting method.}

\item{verbose}{Logical indicating if diagnostic information should be printed.}
}
\description{
Fit a Causal Additive Model (CAM) based on a given DAG
Fit a Causal Additive Model (CAM) based on a given DAG \code{causalDAG} using a specified model 
class (\code{nodeModelName}, \code{nodeModelPars}) for the nodes.
}
\seealso{
fitNode
}
