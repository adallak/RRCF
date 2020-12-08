#source("relaxed_chol.R")
#require(varband)
######################################3
## calculate the correlations on off diagonals
cov.cor<-function(V){
  for(i in 2: nrow(V)){
    for(j in 1:(i-1))
      V[j,i] <- V[i,j] <- V[i,j]/sqrt(V[i,i]*V[j,j])
  }
  diag(V) = 1
  return(V)
}

#################################################
## Calculate penalty of MCP

compareExtended <- function (gl, gt)
{
  ml <- wgtMatrix(ugraph(gl))
  mt <- wgtMatrix(ugraph(gt))
  p <- dim(ml)[2]
  nz_true = sum(mt != 0)
  nz_est  = sum(ml != 0)
  mt[mt != 0] <- rep(1, sum(mt != 0))
  ml[ml != 0] <- rep(1, sum(ml != 0))
  diffm <- ml - mt
  FP = sum(diffm > 0)/2
  nmbTrueGaps <- (sum(mt == 0) - p)/2
  fpr <- if (nmbTrueGaps == 0)
    1
  else (sum(diffm > 0)/2)/nmbTrueGaps
  diffm2 <- mt - ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  TP = nmbTrueEdges - (sum(diffm2 > 0)/2)
  tpr <- if (nmbTrueEdges == 0)
    0
  else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
  trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
  tdr <- if (sum(ml == 1) == 0) {
    if (trueEstEdges == 0)
      1
    else 0
  }
  else trueEstEdges/(sum(ml == 1)/2)
  c(tpr = tpr, fpr = fpr, tdr = tdr, nz_true = nz_true, nz_est = nz_est)
}

#######################
compareL <- function (L_est, L_true)
{
  ml <- L_est
  mt <- L_true
  diag(ml) = diag(mt) = 0
  p <- dim(ml)[2]
  mt[mt != 0] <- rep(1, sum(mt != 0))
  ml[ml != 0] <- rep(1, sum(ml != 0))
  diffm <- ml - mt
  nmbTrueGaps <- (sum(mt == 0) - (p * (p + 1) / 2))
  fpr <- if (nmbTrueGaps == 0)
    1
  else (sum(diffm > 0)) / nmbTrueGaps
  diffm2 <- mt - ml
  nmbTrueEdges <- (sum(mt == 1))
  tpr <- if (nmbTrueEdges == 0)
    0
  else 1 - (sum(diffm2 > 0))/nmbTrueEdges
  trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0))
  tdr <- if (sum(ml == 1) == 0) {
    if (trueEstEdges == 0)
      1
    else 0
  }
  else trueEstEdges/(sum(ml == 1))
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}


###### Function to project on doubly stochastic matrices

proj.ds <- function(P_0, x, y, Z, n.iter = 1000, eps = 1e-8)
{
  k = 1
  n = dim(P_0)[1]
  one.vec = matrix(1, nrow = n)
  converged = FALSE
 # x = y = matrix(0, nrow = n)
#  g = matrix(1 : n , nrow = n)
#  Z = matrix(1, nrow = n, ncol = n)
#  P = P_0 - tcrossprod(x, one.vec) - tcrossprod(one.vec, y) + Z
  ## Initial value for dual variables
  err = c()
  while ((converged == FALSE)&& (k <= n.iter))
  {
    ## Update dual variables
    Z = pmax(matrix(0,n,n), (tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - P_0))
    x = (1 / n) * (P_0 %*% one.vec - (sum(y) + 1) * one.vec + Z %*% one.vec)
    y = (1 / n) * (crossprod(P_0, one.vec) - (sum(x) + 1) * one.vec + crossprod(Z, one.vec))
 #   z = 1/ sum(g^2) * pmax()
    P = P_0 - tcrossprod(x, one.vec) - tcrossprod(one.vec, y) + Z ## Recunstruct primal from dual
    err[k] = abs(primal.obj(P, P_0) - dual.obj(x, y, Z, P_0))
    # cat("dual obj", dual.obj(x, y, Z, P_0), "\n")
    # cat("err obj", k, " and ", err[k], "\n")
    if (err[k] <= eps)
    {
      converged = TRUE
   #   cat("block alogirthm converged", "\n")
    }else{
      k = k + 1
    }
#    cat(" k = ", k, "\n")
#    cat("error is ", err, "\n")
  }
  return = list(P = P, x = x, Z = Z, y = y, err = err )
  return(return)
}

perm.proj <- function(D, rep = 100, S, L, mu, P_old, maj.rule = TRUE)
{
  p = dim(D)[1]
  obj.store = c()
  min = 1e+8
  P.new = P_old # diag(p)
  i = 1
  if(isTRUE(maj.rule))
  {
    while((i <= rep))
  {
    x = runif(p)
    rnk_x = rank(x)
    rnk_Dx = rank(D %*% x)
    ind = match(rnk_Dx, rnk_x)
    P = matrix(0, p, p)
    P[cbind(c(1 : p), ind)] = 1
    store = matrix(0, rep,p)
    obj.store[i] = subobj.function(P, S, L, mu)
    store[i,] = ind
    P.new = P
    i = i + 1
  }
    maj.rule = apply(store, 2, function(x) tabulate(x)/rep )
    P = matrix(0, p, p)
    P[cbind(c(1 : p), maj.rule)] = 1
  }else{
    while((i < rep))
    {
      old.obj = subobj.function(P_old, S, L, mu)
      P = matrix(0, p, p)
      x = runif(p)
      rnk_x = rank(x)
      rnk_Dx = rank(D %*% x)
      ind = match(rnk_Dx, rnk_x)
      P[cbind(c(1 : p), ind)] = 1
      obj.store[i] = subobj.function(P, S, L, mu)
      if  ((obj.store[i] < old.obj))
      {
        min = obj.store[i]
        P.new = P
      }
      i = i + 1
    }
  }
  return(P = P)
}

primal.obj <- function (P, P_0)
{
  obj = 1/2 * sum((P - P_0)^2)
  return(obj)
}

dual.obj <- function(x, y, Z, P_0)
{
  n = length(x)
  one.vec = matrix(1, nrow = n)
  obj = -1 / 2 * sum((tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - Z)^2) - sum(diag(crossprod(Z, P_0))) + crossprod(x, (P_0 %*% one.vec - one.vec)) + crossprod(y, (crossprod(P_0, one.vec) - one.vec))
  return(obj)
}

obj.function <- function(P, S, L, mu, gamma = 2, lambda = 0, penalty = c("lasso", "MCP"))
{
  n <- dim(P)[1]
  one.vec <- matrix(1, n)
  Pii <-  diag(1, n)  - 1/n * tcrossprod(one.vec, one.vec)
#  XPL <- X %*% P %*% L
  if(penalty == "lasso")
  {
    obj <- 1 / 2 * sum(diag(P %*%tcrossprod(S, P) %*% crossprod(L))) - sum(log(diag(L)))  + lambda * sum(abs(L - diag(diag(L))))  # - mu / 2 * sum((Pii %*% P)^2)
     }else{
    row_sum = 0
    for ( i in 2:p)
    {
      index = 1 : i
      row_sum = row_sum + rho(L[1, index], lambda = lambda, gamma = gamma)
    }
    obj <- 1 / 2 * sum(diag(P %*%tcrossprod(S, P) %*% crossprod(L))) - sum(log(diag(L)))  + row_sum # - mu / 2 * sum((Pii %*% P)^2)

    }
  return(obj)
}

subobj.function <- function(P, S, L, mu)
{
  n <- dim(P)[1]
  one.vec <- matrix(1, n)
  Pii <-  diag(1, n)  - 1/n * tcrossprod(one.vec, one.vec)
  #  XPL <- X %*% P %*% L
  obj <- 1 / 2 * sum(diag(P %*%tcrossprod(S, P) %*% crossprod(L)))#  - mu / 2 * sum((Pii %*% P)^2)
  return(obj)
}

gradient <- function(S, P, L, mu)
{
  n <- dim(P)[1]
  one.vec <- matrix(1,nrow = n, ncol = 1)
  Pii <-  diag(1, n) - 1 / n * tcrossprod(one.vec, one.vec)
  grad_f = crossprod(L) %*% P %*% S - mu * Pii %*% P
  return(grad_f)
}


#' This is the main function. Implements the two step RRCF algorithm.
#'
#' @param X - n x p data matrix
#' @param mu - learning rate for the projected gradient algorith
#' @param alpha - used in project gradient algorithm
#' @param s     - learning rate in gradient step
#' @param lambda - penalty hyper-parameter
#' @param gamma  - penalty hyper-parameter. Active when penalty = MCP
#' @param n.iter - global number of itereation
#' @param n.iter.proj - local number of iteration
#' @param eps    - tolerance
#' @param P.init - initial permutation matrix
#' @param penalty - penalty type: MCP or lasso
#'
#' @return permutation matrix P, Cholesky factor L
#' @export
#'
dagrrcf <- function(X, mu = 0.03 , alpha = 1.5, s = 0.01, lambda = 0, gamma = 2, n.iter = 100, n.iter.proj = 100,
                    eps = 1e-4, P.init = NULL,  penalty = c("lasso", "MCP"))
{
  n <- dim(X)[2]
  obs <-dim(X)[1]
  P_proj <- matrix(0, n, n)
#  L <- diag(1, n)
  if (is.null(P.init))
  {
  P <- matrix(rnorm(n^2), n, n)
  }else{
    P <- P.init
  }
  ## initial values for dual variables
  x = y = matrix(1, nrow = n)
  Z = matrix(1, nrow = n, ncol = n)

  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / obs ## Covariance matrix
  L = init = diag(1 / sqrt(diag(S)))   ## Initial value for CSCS
  iter = 1
  converged = FALSE
  err = c()
  while ( (converged == FALSE) && (iter <= n.iter) )
  {
    P_old = P
    L_old = L
    ### Projection into doubly stochastic space
  ## Update P
    P_0 = P - s * gradient(S, P, L, mu)
  ## Projection step
    proj = proj.ds(P_0, x, y, Z, n.iter = n.iter.proj, eps = eps)
    D_proj = proj$P
    ## For warm start
    x = proj$x
    y = proj$y
    Z = proj$Z
    proj_err = proj$err
    P_proj = P_proj + alpha * (D_proj - P_proj) ## Update P

    ## Project into permutation matrix space
    P = perm.proj(P_proj, rep = 50, S, L, mu, P_old)$P

  ### Update L
#  XP_t = tcrossprod(X, P)
  if (penalty == "lasso")
  {
    cv = relaxchol(X = X, P = P, gamma = 0, lambda = lambda,penalty = "lasso", eps = eps)
  }else{
    cv = relaxchol(X = X, P = P, gamma = gamma, lambda = lambda, penalty = "MCP", eps = eps)
  }
  L = cv
  ### Objective update
  interm.obj = obj.function(P, S, L_old, mu, lambda = lambda, gamma = gamma, penalty = penalty)
  old.obj = obj.function(P_old, S, L_old, mu, lambda = lambda, gamma = gamma, penalty = penalty)
  new.obj = obj.function(P, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
  err[iter] = abs(old.obj - new.obj)
  cat("iter error ", err[iter], "\n")
  if (err[iter] < eps)
  {
    converged = TRUE
  }else{
    iter = iter + 1
  }
  }
  return = list(P = P, P_0 = D_proj, L = L, err = err, proj.err = proj_err)
  return (return)
}

chol_lower <- function(omega)
{
  p = dim(omega)[1]
  L = matrix(0, p, p)
  L[p, p] = sqrt(omega[p, p])
  L[p, 1: (p - 1)] = omega[p, 1: (p - 1)] / L[p, p]
  for (i in (p - 1) : (1))
  {
    L[i, i] = sqrt(omega[i, i] - sum(L[c((i + 1) : p), i]^2))
    for ( j in (i - 1) : 1)
    {
      L[i, j] = 1 / L[i, i] * (omega[i, j] - sum(L[c((i + 1): p), i] * L[c((i + 1): p), j]))
    }
  }
  return(L)
}
