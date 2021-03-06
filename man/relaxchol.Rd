% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relaxed_chol.R
\name{relaxchol}
\alias{relaxchol}
\title{Updates the Cholesky factor L, when the permutation matrix is given}
\usage{
relaxchol(X, P, gamma, lambda, maxiter = 100, eps = 1e-04,
  initL = NULL, penalty = c("lasso", "MCP"))
}
\arguments{
\item{X}{- n x p data matrix}

\item{P}{- permutation matrix}

\item{gamma}{- Hyperparameter for penalt. Active when penalty is  MCP.}

\item{lambda}{- Penalty function}

\item{maxiter}{- Maimum number of iterations}

\item{eps}{- tolerance}

\item{initL}{- initial Cholesky Factor}

\item{penalty}{- penaty type: MCP or lasso}
}
\value{
- Cholesky Factor L
}
\description{
Updates the Cholesky factor L, when the permutation matrix is given
}
