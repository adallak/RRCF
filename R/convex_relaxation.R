source("relaxed_chol.R")
require(Matrix)
require(varband)
require(ppcor)
require(graph)
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


###### Function to preject on doubly stochastic matrices


grad.proj<- function(L, P.init = NULL, S,  n.iter = 1000, alpha = 1, s = 0.01, eps = 1e-4,mu = 0.03)
{
  p = dim(S)[1]
  if (is.null(P.init))
  {
    P_proj <- matrix(1/p, p,p)
  }else{
    P_proj <- diag(p)[P.init,]
  }
  ## initial values for dual variables
  x = y = matrix(1, nrow = p)
  Z = matrix(1, nrow = p, ncol = p)
  iter = 1
  converged = FALSE
  err = c()

  while ( (converged == FALSE) && (iter <= n.iter) )
  {
    P_old = P_proj
    L_old = L
    ### Projection into doubly stochastic space
    ## Update P
    P_0 = P_old - s * gradient(S, P_old, L, mu)
    ## Projection step
    proj = proj.ds(P_0, x, y, Z, n.iter = n.iter, eps = eps)
    D_proj = proj$P
    ## For warm start
    x = proj$x
    y = proj$y
    Z = proj$Z
    proj_err = proj$err
    P_proj = P_proj + alpha * (D_proj - P_proj) ## Update P
    err[iter] = sqrt(sum((P_proj - P_old)^2 ))
    if (err[iter] < eps)
    {
      converged = TRUE
    }else{
      iter = iter + 1
    }
  }
  return(list(P = P_proj,err = err))
}


##### Gradient project with lasso penalty ################

sgrad.proj <- function(S, L, p.lambda, s= 0.3, alpha = 1, pmax_iter = 100, eps = 1e-4, P.init = NULL)#, thrsh = TRUE, thrsh.value = 1e-6)
{
  err = c()
  p = dim(S)[1]
  Omega = crossprod(L)
  if (is.null(P.init))
  {
    P =  matrix(1/p, p ,p)
  }else
  {
    P = diag(p)[P.init,]
  }
  x = y = matrix(1, nrow = p)
  Z = matrix(1, nrow = p, ncol = p)
  iter = 1
  converge = FALSE
  P.new = P
  old.obj = subgrad.obj(P, S, L, p.lambda)
  while(converge == FALSE && iter <= pmax_iter)
  {
    P_old = P
    ## Update P
    P_0 = P - s * subgradient(S, P , Omega, p.lambda)
    ## Projection step
    proj = proj.ds(P_0, x, y, Z, n.iter = pmax_iter, eps = eps)
    P.proj = proj$P
    # cat("P here is ", P)
    ## For warm start
    x = proj$x
    y = proj$y
    Z = proj$Z
    proj_err = proj$err
    P.new = P.new + alpha * (P.proj - P.new) ## Update P
    new.obj = subgrad.obj(P.new, S, L, p.lambda)
    if( new.obj < old.obj)
    {
	 P = P.new
   	 err[iter] = sum((old.obj - new.obj)^2)^(1/2)
 	 old.obj = new.obj
   	 if (err[iter] < eps)
    	{
      	converge = TRUE
    	}
    }
   iter = iter + 1
  }
  if(iter == pmax_iter)
  {
    cat("Algorithm does not converged")
  }
 # if (isTRUE(thrsh))
 # {
 #   P[P <= thrsh.value] = 0
 # }
  return(list(P = P, err = err))
}


proj.ds <- function(P_0, x, y, Z, n.iter = 1000, eps = 1e-12)
{
  k = 1
  n = dim(P_0)[1]
  one.vec = matrix(1, nrow = n)
  converged = FALSE
  P_0 = as.matrix(P_0)
  ## Initial value for dual variables
  err = c()
  while ((converged == FALSE)&& (k <= n.iter))
  {
    ## Update dual variables
    Z = pmax(matrix(0,n,n), (tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - P_0))
    x = (1 / n) * (P_0 %*% one.vec - (sum(y) + 1) * one.vec + Z %*% one.vec)
    y = (1 / n) * (crossprod(P_0, one.vec) - (sum(x) + 1) * one.vec + crossprod(Z, one.vec))
    P = P_0 - tcrossprod(x, one.vec) - tcrossprod(one.vec, y) + Z ## Recunstruct primal from dual
    err[k] = abs(primal.obj(P, P_0) - dual.obj(x, y, Z, P_0))
    if (err[k] <= eps)
    {
      converged = TRUE
    }else{
      k = k + 1
    }
  }
  return = list(P = P, x = x, Z = Z, y = y, err = err )
  return(return)
}


sa.proj <- function(D, rep = 100, S , L, mu, P)
{
	T_list = 1 / (1:rep)
	for ( i in 1 : rep)
	{
		P_old = P
		old_obj = subobj.function(P_old, S, L, mu)
		Temp = T_list[i]
    		x = runif(p)
    		rnk_x = rank(x, ties.method = "random")
    		rnk_Dx = rank(D %*% x, ties.method = "random")
    		P = match(rnk_Dx, rnk_x)
      		obj.store = subobj.function(P = P, S, L, mu)
		alpha = min(1, exp(-1 / Temp * (obj.store - old_obj)))
		change = rbinom(1, 1, alpha)
		P = P * change + (1 - change) * P_old
	}
	return(P)
}

perm.proj <- function(D, rep = 100, S, L, mu, P_old, maj.rule = FALSE)
{
  p = dim(D)[1]
  obj.store = c()
  i = 1
 # old.obj = subobj.function(P_old, S, L, mu)
  if(isTRUE(maj.rule))
  {
    store = matrix(0, rep, p)
    store[rep,] = p
    while((i < rep))
  {
    x = runif(p)
    rnk_x = rank(x, ties.method = "random")
    rnk_Dx = rank(D %*% x, ties.method = "random")
    ind = match(rnk_Dx, rnk_x)
    store[i,] = ind
    i = i + 1
    }
    h_x = apply(store, 2, function(x) tabulate(x))

    E_x = rep / p
    s = 2^(- E_x/ h_x)
    temp = apply(s, 2, function(x) which.max(x))
#    cat("Initial temp is", temp, "\n")
    k = 1 : p
    not_included = k[-which(1:p %in% temp)]
 #   cat("not included is", not_included, "\n")
    dupl_score = s[,which(duplicated(temp))]
    duplic_index =  which(duplicated(temp))
    for(i in not_included)
    {
      if(!is.null(ncol(dupl_score)))
      {
        l = which.max(dupl_score[i,])
        temp[duplic_index][l] = i
        dupl_score = dupl_score[, -l]
        duplic_index = duplic_index[-l]
      }else{
        temp[duplic_index] = i
      }
    }
    rm(store,not_included, duplic_index, dupl_score, h_x, s)
    gc()
    P = temp
  }else{

    #old.obj = subobj.function(P_old, S, L, mu)
    min = 1e+8
    P = P_old
    while((i < rep))
    {
      x = runif(p)
      rnk_x = rank(x, ties.method = "random")
      rnk_Dx = rank(D %*% x, ties.method = "random")
      ind = match(rnk_Dx, rnk_x)
      obj.store = subobj.function(P = ind, S, L, mu)
      if  ((obj.store < min))
      {
        min = obj.store
        P = ind
      }
      i = i + 1
    }
  }
  return(P)
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
  obj = -1 / 2 * sum((tcrossprod(x, one.vec) + tcrossprod(one.vec, y) - Z)^2) - sum(diag(crossprod(Z, P_0))) +
    crossprod(x, (P_0 %*% one.vec - one.vec)) + crossprod(y, (crossprod(P_0, one.vec) - one.vec))
  return(obj)
}

obj.function <- function(P, S, L, mu, gamma = 2, lambda = 0, penalty = c("lasso", "MCP"))
{
  p <- dim(S)[1]
  one.vec <- matrix(1, p)
  Pii <-  diag(1, p)  - 1/p * tcrossprod(one.vec, one.vec)
  if(penalty == "lasso")
  {
    obj <- 1 / 2 * sum(diag(S[P,P] %*% crossprod(L))) - sum(log(diag(L)))  + lambda * sum(abs(L - diag(diag(L))))  # - mu / 2 * sum((Pii %*% P)^2)
     }else{
    row_sum = 0
    for ( i in 2:p)
    {
      index = 1 : i
      row_sum = row_sum + rho(L[1, index], lambda = lambda, gamma = gamma)
    }
    obj <- 1 / 2 * sum(diag(S[P,P] %*% crossprod(L))) - sum(log(diag(L)))  + row_sum # - mu / 2 * sum((Pii %*% P)^2)

    }
  return(obj)
}

subobj.function <- function(P, S, L, mu)
{
  p <- dim(S)[1]
  one.vec <- matrix(1, p)
  Pii <-  diag(1, p)  - 1/p * tcrossprod(one.vec, one.vec)
  obj <- 1 / 2 * sum(diag(S[P,P] %*% crossprod(L)))#  - mu / 2 * sum((Pii %*% P)^2)
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

subgradient <- function(S, P , Omega, lambda)
{
  subgrad = sign(P)
 # subgrad[subgrad == 0] = runif(length(subgrad[subgrad == 0]), -1, 1)
  return(Omega %*% P %*% S + lambda * subgrad)
}

subgrad.obj <- function(P, S, L, lambda)
{
  obj <- 1 / 2 * sum(diag(P %*% S %*% t(P) %*% crossprod(L))) + lambda * sum(abs(P))  #  - mu / 2 * sum((Pii %*% P)^2)
  return(obj)
}
# Armijo_rule <- function(P, proj_P, beta = 1/5, sigma = 10^(-2), S, L, mu, lambda){
#   d = proj_P - P
#   m = 1
#   alpha = beta
#   conv = obj.function(P, S, L, mu, lambda) - obj.function(d, S, L, mu, lambda)
#   #cat("conv is ", conv,"\n")
#   while (conv <= - sigma * beta^m * sum(diag(crossprod(gradient(S, P, L, mu),d)))){
#     alpha = beta^m
#     #cat("difference is ", conv  - (- sigma * beta^m * sum(diag(crossprod(gradient(S, P, L, mu),d)))), "\n")
#     m = m + 1
#   }
#   return(alpha)
# }

dag.cr <- function(X, mu = 0.03 , alpha = 1, s = 0.01, lambda = 0, gamma = 2, n.iter = 100, n.iter.proj = 100,
                   eps = 1e-4, penalty = c("lasso", "MCP"), maj.rule = FALSE, perm.rep = 100, refine = FALSE, ref.alpha = 0.001, BH = TRUE)
{
  p <- dim(X)[2]
  obs <-dim(X)[1]
  P_proj <- matrix(0, p, p)
#  L <- diag(1, n)
#  if (is.null(P.init))
#  {
 # 	P <- 1:p
#  }else{
 #   	P <- P.init
#  }
  ## initial values for dual variables
  x = y = matrix(1, nrow = p)
  Z = matrix(1, nrow = p, ncol = p)
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
    ord = perm.proj(P_proj, rep = perm.rep, S, L, mu, P, maj.rule = maj.rule)
    P = diag(p)[ord, ]
    cat("order is \n")
    print(ord)
  ### Update L
  if (penalty == "lasso")
  {
    L = relaxchol(X = X, P = ord, gamma = 0, lambda = lambda,penalty = "lasso", eps = eps)
  }else{
    L = relaxchol(X = X, P = ord, gamma = gamma, lambda = lambda, penalty = "MCP", eps = eps)
  }
  ### Objective update
#  interm.obj = obj.function(P, S, L_old, mu, lambda = lambda, gamma = gamma, penalty = penalty)
  new.obj = obj.function(P, S, L, mu, lambda = lambda, gamma = gamma, penalty = penalty)
  err[iter] = abs(old.obj - new.obj)
  cat("iter error ", err[iter], "\n")
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
  return = list(P = P, P_0 = D_proj, B = B,ord = ord, B.p = B.p, L = L, err = err, proj.err = proj_err)
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


# dag.cr_CV <- function(X, mu = 0.03 , alpha = 1, s = 0.01, lamlist = NULL, nlam = 60,  flmin = 1e-2, folds = NULL, nfolds = 5, gamma = 2, n.iter = 100, n.iter.proj = 100, eps = 1e-4, P.init = NULL,  penalty = c("lasso", "MCP"))
# {
#   n <- dim(X)[2]
#   obs <-dim(X)[1]
#   #  L <- diag(1, n)
#   if (is.null(P.init))
#   {
#     P <- matrix(rnorm(n^2), n, n)
#   }else{
#     P <- P.init
#   }
#   ## initial values for dual variables
#   x = y = matrix(0, nrow = n)
#   Z = matrix(0, nrow = n, ncol = n)
#
#   S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / obs ## Covariance matrix
#   L = init = diag(1 / sqrt(diag(S)))   ## Initial value for CSCS
#  # old.obj = obj.function(P, S, init, mu, lambda = lambda)
# #  new.obj = 100
#   converged = FALSE
#   iter = 0
#   err = c()
# #  while ( (iter <= n.iter) )
#   while((converged == FALSE) && (iter < n.iter))
#    {
#     P_old = P
#   #  new.obj = old.obj
#     ## Update P
#     P_0 = P - s * gradient(S, P, L, mu)
#     ## Projection step
#     proj = proj.ds(P_0, x, y, Z, n.iter = n.iter.proj, eps = eps)
#     D_proj = proj$P
#     P_proj = perm.proj(D_proj, rep = 1000, S, L, mu)
#     # cat("P_proj", P_proj)
#     ## For warm start
#     x = proj$x
#     y = proj$y
#     Z = proj$Z
#     proj_err = proj$err
#     P = P + alpha * (P_proj - P) ## Update P
#
#     ### Update L
# #    XP_t = tcrossprod(X, P)
#     if (penalty == "lasso")
#     {
#       cv = relaxchol_CV(X = X , gamma = 0, P = P, lamlist = lamlist, nlam = nlam, flmin = flmin, folds = folds, nfolds = nfolds, maxiter = n.iter, penalty = c("lasso"), eps = eps)
#
#     }else{
#       cv = relaxchol_CV(X = X , gamma = gamma, P = P, lamlist = lamlist, nlam = nlam, flmin = flmin, folds = folds, nfolds = nfolds, maxiter = n.iter, penalty = c("MCP"), eps = eps)
#     }
#     #  S_perm <- P %*% tcrossprod(S, P)
#     #S_perm = P %*% tcrossprod(S, P) ## Permute the matrix
#     #  cat("S is", S_perm, "\n")
#     #  L = varband(S_perm, lambda = lambda, init, w = FALSE, lasso = TRUE)
#     #  init <- L
#     L = cv$L_fit
#     lambda = cv$lambda_cv
#     #lambda = cv$lamlist[cv$ibest_fit]
#     #}
#     ### Objective update
#  #   cat("iter error ", err[iter], "\n" )
#     if(abs(obj.function(P, S, L, mu, lambda = lambda) - obj.function(P_old, S, L, mu, lambda = lambda))< eps)
#     {
#       converged = TRUE
#     }else{
#     iter = iter + 1
#     }
#     # if (err[iter] < eps)
#     # {
#     #   print("Iteration Converges")
#     #   break
#     # }
#   }
#   return = list(P = P, P_0 = P_0, L = L, lambda = lambda, proj.err = proj_err)
#   return (return)
# }




