% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convex_relaxation.R
\name{dagrrcf}
\alias{dagrrcf}
\title{This is the main function. Implements the two step RRCF algorithm.}
\usage{
dagrrcf(X, mu = 0.03, alpha = 1, s = 0.01, lambda = 0, gamma = 2,
  n.iter = 100, n.iter.proj = 100, eps = 1e-04,
  penalty = c("lasso", "MCP"), perm.rep = 100, refine = FALSE,
  ref.alpha = 0.001, BH = TRUE)
}
\arguments{
\item{X}{- n x p data matrix}

\item{mu}{- learning rate for the projected gradient algorith}

\item{alpha}{- used in project gradient algorithm}

\item{s}{- learning rate in gradient step}

\item{lambda}{- penalty hyper-parameter}

\item{gamma}{- penalty hyper-parameter. Active when penalty = MCP}

\item{n.iter}{- global number of itereation}

\item{n.iter.proj}{- local number of iteration}

\item{eps}{- tolerance}

\item{penalty}{- penalty type: MCP or lasso}

\item{P.init}{- initial permutation matrix}
}
\value{
permutation matrix P, Cholesky factor L
}
\description{
This is the main function. Implements the two step RRCF algorithm.
}
