# RRCF
## Introduction

The `RRCF` is a package that learns Bayesian Networks from the Birkhoff polytope. 

This document serves as an introduction of using the package.

The main function is `dagrrcf`, which takes a data matrix of the observations and hyper-parameters to return the permutation matrix and Cholesky factor. 

## Installation

To install the latest version from Github, use

```s
library(devtools)
devtools::install_github("adallak/RRCF")
```
Simple example

```s
require(RRCF)
p = 10
n = 150
prob = 0.3
X = genDAGdata(n, p, prob)
est_P = dagrrcf(X, mu = 0.1 , alpha = 1, s = 2, lambda = 0.1, gamma = 2, n.iter = 100, penalty = c("MCP"), perm.rep = 100)$P
```
