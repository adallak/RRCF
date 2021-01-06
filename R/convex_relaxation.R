#source("relaxed_chol.R")
#require(Matrix)
#require(varband)
#require(ppcor)
#require(graph)
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
  c(tpr = tpr, fpr = fpr, tdr = tdr, nz_est = nz_est)
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


# proj.ds <- function(P_0, x, y, Z, n.iter = 1000, eps = 1e-12)
# {
#   k = 1
#   n = dim(P_0)[1]
#   one.vec = matrix(1, nrow = n)
#   converged = FALSE
#   P_0 = as.matrix(P_0)
#   ## Initial value for dual variables
#   err = c()
#   while ((converged == FALSE)&& (k <= n.iter))
#   {
#     ## Update dual variables
#     Z = pmax(matrix(0,n,n), (tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - P_0))
#     x = (1 / n) * (P_0 %*% one.vec - (sum(y) + 1) * one.vec + Z %*% one.vec)
#     y = (1 / n) * (crossprod(P_0, one.vec) - (sum(x) + 1) * one.vec + crossprod(Z, one.vec))
#     P = P_0 - tcrossprod(x, one.vec) - tcrossprod(one.vec, y) + Z ## Recunstruct primal from dual
#     err[k] = abs(primal.obj(P, P_0) - dual.obj(x, y, Z, P_0))
#     if (err[k] <= eps)
#     {
#       converged = TRUE
#     }else{
#       k = k + 1
#     }
#   }
#   return = list(P = P, x = x, Z = Z, y = y, err = err )
#   return(return)
# }

proj.ds <- function(P_0, x, y, Z, z, n.iter = 1000, eps = 1e-12)
{
      k = 1
      n = dim(P_0)[1]
      one.vec = matrix(1, nrow = n)
      converged = FALSE
      P_0 = as.matrix(P_0)
      ## Initial value for dual variables
      err = c()

      D =   matrix(0, n, 1)
      D[1] = 1
      D[n] = -1
      delta = 1
      g = matrix(1 : n, n , 1)

      while ((converged == FALSE)&& (k <= n.iter))
      {
            ## Update dual variables
            Z = pmax(matrix(0,n,n), (tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - P_0 + z * D %*% t(g)))
            x = (1 / n) * (P_0 %*% one.vec - (sum(y) + 1) * one.vec + Z %*% one.vec - z * D * as.numeric(crossprod(g, one.vec)))
            y = (1 / n) * (t(P_0) %*% one.vec - (sum(x) + 1) * one.vec +t(Z) %*% one.vec - z * g * as.numeric(crossprod(D, one.vec)))
            z = (1 / sum(g^2)) * max(0, (crossprod(D, (Z + P_0)) %*% g + delta - crossprod(D, x) * as.numeric(crossprod(g, one.vec)) +
                                               crossprod(D, one.vec) * as.numeric(crossprod(g, y))) / crossprod(D) )
            P = P_0 - tcrossprod(x, one.vec) - tcrossprod(one.vec, y) + Z  - z * D %*% t(g)## Recunstruct primal from dual
            err[k] = abs(primal.obj(P, P_0) - dual.obj(x, y, Z, P_0, z))
            #          cat("inner error is ", err[k], "\n")
            if (err[k] <= eps)
            {
                  converged = TRUE
            }else{
                  k = k + 1
            }
      }
      return = list(P = P, x = x, Z = Z, y = y, err = err, z = z )
      return(return)
}


perm.proj <- function(D, rep = 100, S, L, mu, P_old)
{
  p = dim(D)[1]
  obj.store = c()
  i = 1
 # old.obj = subobj.function(P_old, S, L, mu)
   #old.obj = subobj.function(P_old, S, L, mu)
  min = 1e+8
  P = P_old
  while(i < rep)
  {
      x = sort(runif(p))
        #            rnk_x = rank(x, ties.method = "random")
      rnk_Dx = rank(D %*% x, ties.method = "random")
      ind = rnk_Dx
      P_current = diag(p)[,ind]
      obj.store = subobj.function(P = P, S, L, mu)
      if  (obj.store < min)
      {
        min = obj.store
        P = P_current
      }
      i = i + 1
   }
  return(P)
}

primal.obj <- function (P, P_0)
{
  obj = 1/2 * sum((P - P_0)^2)
  return(obj)
}

# dual.obj <- function(x, y, Z, P_0)
# {
#   n = length(x)
#   one.vec = matrix(1, nrow = n)
#   obj = -1 / 2 * sum((tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - Z)^2) - sum(diag(crossprod(Z, P_0))) +
#     crossprod(x, (P_0 %*% one.vec - one.vec)) + crossprod(y, (crossprod(P_0, one.vec) - one.vec))
#   return(obj)
# }
dual.obj <- function(x, y, Z, P_0 , z)
{
      n = length(x)
      one.vec = matrix(1, nrow = n)
      D  =  matrix(0, n, 1)
      D[1] = 1
      D[n] = -1
      delta = 1
      g = matrix(1 : n, n ,1 )
      obj = -1 / 2 * sum((tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - Z + z * D %*% t(g))^2) - sum(diag(crossprod(Z, P_0))) +
            crossprod(x, (P_0 %*% one.vec - one.vec)) + crossprod(y, (crossprod(P_0, one.vec) - one.vec)) +
            z * (crossprod(D, P_0) %*% g + delta)
      return(as.numeric(obj))
}

obj.function <- function(P, S, L, mu, gamma = 2, lambda = 0, penalty = c("lasso", "MCP"))
{
      p <- dim(S)[1]
      one.vec <- matrix(1, p)
      Pii <-  diag(1, p)  - 1/p * tcrossprod(one.vec, one.vec)
      if(penalty == "lasso")
      {
            obj <- 1 / 2 * sum(diag(P %*% S %*% t(P) %*% crossprod(L))) - sum(log(diag(L)))  +
                  lambda * sum(abs(L - diag(diag(L)))) - mu / 2 * sum((Pii %*% P)^2)
      }else{
            row_sum = 0
            for ( i in 2:p)
            {
                  index = 1 : i
                  row_sum = row_sum + rho(L[1, index], lambda = lambda, gamma = gamma)
            }
            obj <- 1 / 2 * sum(diag(P %*% S %*% t(P) %*% crossprod(L))) - sum(log(diag(L)))  +
                  row_sum  - mu / 2 * sum((Pii %*% P)^2)

      }
      return(obj)
}

subobj.function <- function(P, S, L, mu)
{
      p <- dim(S)[1]
      one.vec <- matrix(1, p)
      Pii <-  diag(1, p)  - 1/p * tcrossprod(one.vec, one.vec)
      obj <- 1 / 2 * sum(diag(P %*% S %*% t(P) %*% crossprod(L))) - mu / 2 * sum((Pii %*% P)^2)
      return(obj)
}

gradient <- function(S, P, L, mu)
{
  p <- dim(S)[1]
  one.vec <- matrix(1,nrow = p, ncol = 1)
#  P.mat = as(as.integer(P),"pMatrix")
#  P = diag(p)[P,]
  Pii <-  diag(1, p) - 1 / p * tcrossprod(one.vec, one.vec)
  grad_f = crossprod(L) %*% P%*% S - mu * Pii %*% P
  return(grad_f)
}

updateP<- function(x, y, Z,z, niter, s, S, L, mu, eps, alpha, lambda, gamma, penalty, perm.rep = 500, P.init = NULL)
{
       if (is.null(P.init))
       {
             P_proj = matrix(1 / p, p, p)
       }else{
             P_proj = P.init
       }
      iter = 1
      converged = FALSE
      err = c()
      while((converged == FALSE) & (iter < niter))
      {
            s = s / iter
            P_old = P_proj
            P_0 = P_old - s * gradient(S, P_old, L, mu)
            ## Projection step
            proj = proj.ds(P_0, x, y, Z, z, n.iter = 5000, eps = eps)
            D_proj = proj$P
            ## For warm start
            x = proj$x
            y = proj$y
            Z = proj$Z
            z = proj$z
            P_proj = P_old + alpha * (D_proj - P_old) ## Update P
            old.obj = obj.function(P_old, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
            new.obj = obj.function(P_proj, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
            #          cat("error is ",abs(old.obj - new.obj), "\n" )
            if (abs(old.obj - new.obj) < eps)
            {
                  converged = TRUE
            }else{
                  iter = iter + 1
            }
      }
      ## Project into permutation matrix space
      P = perm.proj(P_proj, rep = perm.rep, S, L, mu, P_old)
      #P = diag(p)[ord, ]
      return(P)
}

# updateP<- function(x, y, Z, niter, s, S, L, mu, eps, alpha, lambda, gamma, penalty, P.init = NULL, perm.rep = 100)
# {
#       if (is.null(P.init))
#       {
#             P_proj = matrix(1 / p, p, p)
#       }else{
#             P_proj = P.init
#       }
#       iter = 1
#       converged = FALSE
#       err = c()
#       while((converged == FALSE) & (iter < niter))
#       {
#             s = s / iter
#             P_old = P_proj
#             P_0 = P_old - s * gradient(S, P_old, L, mu)
#             ## Projection step
#             proj = proj.ds(P_0, x, y, Z, n.iter = 1000, eps = eps)
#             D_proj = proj$P
#             ## For warm start
#             x = proj$x
#             y = proj$y
#             Z = proj$Z
#             P_proj = P_proj + alpha * (D_proj - P_proj) ## Update P
#             old.obj = obj.function(P_old, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
#             new.obj = obj.function(P_proj, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
#             if (abs(old.obj - new.obj) < eps)
#             {
#                   converged = TRUE
#             }else{
#                   iter = iter + 1
#             }
#       }
#       ## Project into permutation matrix space
#       ord = perm.proj(P_proj, rep = perm.rep, S, L, mu, P)
#       P = diag(p)[ord, ]
#       return(list(P = P, ord = ord))
# }


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
#' @example
#' p = 10
#' n = 150
#' prob = 0.3
#' X = genDAGdata(n, p, prob)
#' est_P = dagrrcf(X, mu = 0.1 , alpha = 1, s = 2, lambda = 0.1, gamma = 2, n.iter = 100, n.iter.proj = 100, eps = 1e-4, penalty = c("lasso"), perm.rep = 100, refine = FALSE,ref.alpha = 0.001, BH = TRUE)$P


dagrrcf <- function(X, mu = 0.03 , alpha = 1, s = 0.01, lambda = 0, gamma = 2, n.iter = 100, n.iter.proj = 100,
                   eps = 1e-4, penalty = c("lasso", "MCP"), perm.rep = 100, refine = FALSE,
                   ref.alpha = 0.001, BH = TRUE)
{
  p <- dim(X)[2]
  obs <-dim(X)[1]
#  P_proj <- matrix(0, p, p)
#  L <- diag(1, n)
#  if (is.null(P.init))
#  {
 # 	P <- 1:p
#  }else{
 #   	P <- P.init
#  }
  ## initial values for dual variables
  x = y = matrix(1, nrow = p, 1)
  Z = matrix(1, nrow = p, ncol = p)
  z = 1
  P = matrix(1 / p, p, p)
  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / obs ## Covariance matrix
  L = init = diag(1 / sqrt(diag(S)))   ## Initial value for CSCS
  iter = 1
  converged = FALSE
  old.obj = 1e+8
  err = c()
  while ( (converged == FALSE) && (iter <= n.iter) )
  {
  ### Projection into doubly stochastic space
  ## Update P
        temp = updateP(x, y, Z, z, n.iter, s, S, L, mu, eps, alpha, lambda, gamma, penalty, perm.rep, P)
        P = temp
#        ord = temp$ord
        ### Update L
       if (penalty == "lasso")
       {
            L = relaxchol(X = X, P = P, gamma = 0, lambda = lambda,penalty = "lasso", eps = eps)
       }else{
            L = relaxchol(X = X, P = P, gamma = gamma, lambda = lambda, penalty = "MCP", eps = eps)
       }
        ### Objective update
#  interm.obj = obj.function(P, S, L_old, mu, lambda = lambda, gamma = gamma, penalty = penalty)
      new.obj = obj.function(P, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
      err[iter] = abs(old.obj - new.obj)
        if (err[iter] < eps)
      {
            converged = TRUE
      }else{
      iter = iter + 1
      P_old = P
      L_old = L
      old.obj = obj.function(P_old, S, L_old, mu, lambda = lambda, gamma = gamma, penalty = penalty)

      }
  }
  B.p = - L / matrix(diag(L), p,p,  byrow = TRUE)
  diag(B.p) = 0
  if(isTRUE(refine))
  {
    B = refine(X, P, L, ref.alpha, BH = BH)
  }else{
    L.p = t(P) %*% L %*% P
    B = - L.p / matrix(diag(L.p), p,p,  byrow = TRUE)
    diag(B) = 0
    rm(L.p)
  }
  return = list(P = P, B = B, B.p = B.p, L = L, err = err)
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

refine <- function(X, P, L, alpha, BH = TRUE)
{
  n = dim(X)[1]
  p = dim(L)[1]
  L.p = t(as(as.integer(c(P)),"pMatrix")) %*% L %*% as(as.integer(c(P)),"pMatrix")
  print("L.p is ",L.p)

  B = - as.matrix(L.p) / matrix(diag(L.p), p, p, byrow = TRUE)
  diag(B) = 0
  Z_alpha = qnorm(1 - alpha/2)
  rm(L.p)
  gc()
  for (j in 1 : p)
  {
    p.val = c()
    pa_ind = which(B[j,] != 0)
    pa_j = B[j, pa_ind]
    if (length(pa_ind) == 0)
    {
      next
    }else if (length(pa_ind) == 1){
      p.val = cor(X[,j],X[,pa_ind])
      }else{
      for (k in 1: length(pa_ind))
      {
        s = pa_ind[-k]
        p.val[k] = pcor.test(X[,j], X[,k], X[,s])$p.value
      }
      if (isTRUE(BH))
      {
        p.val = p.adjust(p.val, method = "BH")
      }
      remov.ind = which(p.val < alpha)
      B[remov.ind, j] = 0
    }
  }
  return(B)
}


########################################3
#### This function are borrowed from https://github.com/WY-Chen/EqVarDAG

randomDAG2_er <- function(p,probConnect)
{
  # This function is modified from randomDAG2 function by Jonas Peters
  DAG <- diag(rep(0,p))
  causalOrder <- sample(p)
  for(i in 3:(p))
  {
    node <- causalOrder[i]
    possibleParents <- causalOrder[1:(i-1)]
    Parents <- possibleParents[rbinom(length(possibleParents),1,probConnect)==1]
    DAG[Parents,node] <- rep(1,length(Parents))
  }
  node <- causalOrder[p-1]
  ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
  DAG[causalOrder[1],causalOrder[2]] <- 1

  return(list(DAG=DAG,TO=causalOrder))
}
# Generate data from SEM
#' Title
#'
#' @param n - number of observations
#' @param p - number of variables
#' @param pc - probability of non zero edges
#'
#' @return n x p data matrix X,  p x p true coefficent matrix B.p
#' @export
#'
#' @examples
#' n = 100
#' p = 10
#' pc = 0.2
#' X = genDAGdata(n, p, pc)
genDAGdata<-function(n,p,pc){
  Bmin<-0.3
  D<-randomDAG2(p,pc)
  truth<-D$DAG
  TO<-D$TO
  errs <- matrix(rnorm(p * n), nrow = p, ncol = n)
  B<-t(truth)
  B[B==1]<-runif(sum(truth),Bmin,1)*(2*rbinom(sum(truth),1,0.5)-1)
  X <- solve(diag(rep(1, p)) - B, errs)
  X <- t(X)
  B.p = B[TO, TO]
  return(list(truth=truth,B=B, B.p = B.p, X=X,TO=TO))
}


