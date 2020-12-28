
soft <- function(a, lambda){
  return(as.numeric(sign(a) * pmax(0, abs(a) - lambda)))
}

rho <- function(x, lambda, gamma)
{
  n = length(x)
  result = 0
  for (i in 1:(n-1))
  {
    if(abs(x[i]) < lambda * gamma)
    {
      result = result + (lambda * abs(x[i]) - (x[i]^2)/ (2 * gamma))
    } else {
      result = result + gamma * (lambda^2) / 2
    }
  }
  return(result)
}

h_x <- function(x, A, lambda, gamma,penalty = c("lasso", "MCP"))
{
  k = length(x)
  if (penalty == "MCP")
  {
    h = -2 * log(x[k]) + crossprod(x, A) %*% x + rho(x, lambda = lambda, gamma = gamma)
  }else{
    h = -2 * log(x[k]) + crossprod(x, A) %*% x + lambda * sum(abs(x[-k]))
  }
  return(h)
}


RC_i <- function(i, A, gamma, lambda, maxiter = 100, eps = 1e-4, initx = NULL, penalty = c("lasso", "MCP"))
{
  if (i == 1)
  {
    x = 1 / sqrt(A)
  }else{
    p = dim(A)[1]
    iter = 1
    converged = FALSE
    if (is.null(initx))
    {
      x = matrix(0, nrow = i, ncol = 1)
      x[i] = sqrt(A[i,i])
    }else{
      x = initx
    }
    while((converged == FALSE) && (iter < maxiter))
    {
      x_old = x
      for (j in 1 : (i-1))
      {
        A_jj = A[j, j]
        y_j = sum(A[-j, j] * x[-j])
        if(penalty == "MCP")
        {
          if ((abs(y_j) / A_jj) <= (gamma * lambda))
          {
            x[j] = soft((-2 * y_j), lambda) / (2 *  A_jj - 1/ gamma)
          }else{
            x[j] = - y_j / A_jj
          }
        }else{
          x[j] = soft((-2 * y_j), lambda) / (2 *  A_jj)
        }
      }
      ### Update Diagonal Element
      A_ii = A[i, i]
        y_i = sum(A[-i, i] * x[-i])
      x[i] = (-y_i + sqrt(y_i^2 + 4 * A_ii)) / (2* A_ii)
      if(max(abs(x - x_old)) < eps)
      {
        converged = TRUE
      }else{
        iter = iter + 1
      }
    }
  }
  return(x)
}

#' Updates the Cholesky factor L, when the permutation matrix is given
#'
#' @param X - n x p data matrix
#' @param P - permutation matrix
#' @param gamma - Hyperparameter for penalt. Active when penalty is  MCP.
#' @param lambda - Penalty function
#' @param maxiter - Maimum number of iterations
#' @param eps     - tolerance
#' @param initL   - initial Cholesky Factor
#' @param penalty - penaty type: MCP or lasso
#'
#' @return  - Cholesky Factor L
#' @export
relaxchol <- function(X, P, gamma, lambda, maxiter = 100, eps = 1e-4, initL = NULL, penalty = c("lasso", "MCP"))
{
  n = dim(X)[1]
  p = dim(X)[2]
  XP = X[,P]# tcrossprod(X, P)
  colnames(XP) = colnames(X)
  S = crossprod(scale(XP, center = TRUE, scale = FALSE)) / n
  if (is.null(initL))
  {
    L = diag(diag(S), p, p)
  }else{
    L = initL
  }
  colnames(L) = rownames(L) = colnames(X)[P]
#  L_cscs = varband(S, lambda = lambda, init = L, lasso = TRUE, w = FALSE )
  for (i in 1 : p)
  {
    index = 1 : i
    A = S[index, index]
    if (penalty == "MCP")
    {
      L[i , index] = RC_i(i, A = A, gamma = gamma, lambda = lambda, maxiter = maxiter, eps = eps, initx = L[i, index], penalty = "MCP")
    }else{
      L[i , index] = RC_i(i, A = A, gamma = 0, lambda = lambda, maxiter = maxiter, eps = eps, initx = L[i, index], penalty = "lasso")
    }
  }
  return(L)
}

##### Tuning using BIC
relaxchol_BIC <- function(X , gamma = 2, P , lamlist = NULL, nlam = 60, mu, flmin = 1e-2, maxiter = 100, penalty = c("lasso", "MCP"), eps = 1e-4, ebic_gam = 0.25)
{
  n <- nrow(X)
  p <- ncol(X)

  S <- crossprod(scale(X, center=TRUE, scale=FALSE)) / n
  initL = diag(diag(S))
  if (is.null(lamlist)) {
    lam_max <- lammax(S = S)
    lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                       flmin = flmin, S = S)
  } else {
    nlam <- length(lamlist)
  }
  theta <- expand.grid(lamlist, gamma)

  bic.fit <- ebic.fit <- c()
  storage <- matrix(0, nrow = p, ncol = p)
  for (i in 1 : dim(theta)[1])
  {
    storage =  relaxchol(X = X, P = P, gamma = theta[i,2], lambda = theta[i,1],
                         initL = initL, maxiter = maxiter, eps = eps, penalty = penalty)
    initL = storage
    bic.fit[i] = n/2 * sum(diag(S[P,P] %*% crossprod(initL))) - n * determinant(initL, logarithm = TRUE)$modulus[1] + sum(initL != 0) * log(max(n,p))
    ebic.fit[i] = n/2 * sum(diag(S[P,P] %*% crossprod(initL))) - n * determinant(initL, logarithm = TRUE)$modulus[1] + sum(initL != 0) * log(n) + 4 * sum(initL != 0) * log(p) * ebic_gam
  }
  i_bic = which.min(bic.fit)
  i_ebic = which.min(ebic.fit)
  gamma.bic = theta[i_bic, 2]
  lambda.bic = theta[i_bic, 1]
  gamma.ebic = theta[i_ebic, 2]
  lambda.ebic = theta[i_ebic, 1]
  return(list(gamma = gamma.bic, lambda = lambda.bic, bic.fit = bic.fit, lambda.ebic = lambda.ebic,gamma.ebic = gamma.ebic ))
}

#### Tuning using CV


relaxchol_path <- function(X, P, gamma, maxiter = 100, eps = 1e-4, penalty = c("lasso", "MCP"),
                         lamlist = NULL, nlam = 60, flmin = 0.01){
  n = dim(X)[1]
  p = dim(X)[2]
  S = crossprod(scale(X, center = TRUE, scale = FALSE)) / n

  if (is.null(lamlist)) {
    lam_max <- lammax(S = S)
    lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                       flmin = flmin, S = S)
  } else {
    nlam <- length(lamlist)
  }
  theta = expand.grid(lamlist, gamma)
  ntheta = dim(theta)[1]
  result<- array(NA, c(p, p, ntheta))

  for (i in seq(ntheta)) {
    if(i==1){
      result[, , i] <- diag(1/sqrt(diag(S)))
    }
    else
   {
      result[, , i] <- relaxchol(X = X, P = P, gamma = theta[i, 2] , lambda = theta[i, 1], maxiter = maxiter, eps = eps, initL = result[, , i - 1], penalty = penalty)

    }
  }
  return(list(path = result, lamlist = lamlist))
}

lammax <- function(S){
  #### This function calculates the max value in the tuning parameter list
  # such that the estimator L_{\lambda} is a diagonal matrix
  # NOTE: this is not necessarily true, but generally
  # a upper bound of the value we are looking for.

  # Args:
  #     S: the p-by-p sample covariance matrix

  p <- ncol(S)
  sighat <- rep(NA, p-1)
  for (r in seq(2, p)){
    sighat[r-1] <- max(abs(S[(1:(r-1)), r]))/sqrt(S[r, r])
  }
  2 * max(sighat)
}

pathGen <- function(nlam, lam_max, flmin, S){
  # Generate a path of lambda, with
  # nlam/2 decreasing exponentially
  # nlam/2 decreasing linearly
  # lam_max <- lammax(S)
  lamlist_lin <- lam_max * exp(seq(0, log(flmin), length = nlam/2))
  lamlist_exp <- seq(lam_max - 1e-8, lam_max*flmin - 1e-8, length.out = nlam/2)
  return(sort(unique(c(lamlist_lin, lamlist_exp)), decreasing = T))
}


relaxchol_CV <- function(X , gamma, P , lamlist = NULL, nlam = 60, flmin = 1e-2, folds = NULL, nfolds = 5, maxiter = 100, penalty = c("lasso", "MCP"), eps = 1e-4)
{
  n <- nrow(X)
  p <- ncol(X)

  S <- crossprod(scale(X, center=TRUE, scale=FALSE)) / n
  if(is.null(folds))
    folds <- makefolds(n, nfolds = nfolds)
  nfolds <- length(folds)

  if (is.null(lamlist)) {
    lam_max <- lammax(S = S)
    lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                     flmin = flmin, S = S)
  } else {
    nlam <- length(lamlist)
  }
  theta = expand.grid(lamlist, gamma)
  ntheta = dim(theta)[1]
  errs_fit <- matrix(NA, ntheta, nfolds)

  # error function is the negative log Gaussian likelihood
  for (i in seq(nfolds)) {
    # train on all but i-th fold:
    x_tr <- X[-folds[[i]],]
    meanx <- colMeans(x_tr)
    x_tr <- scale(x_tr, center = meanx, scale = FALSE)
    S_tr <- crossprod(x_tr) / (dim(x_tr)[1])

    path_fit <- relaxchol_path(X = x_tr, P = P, gamma = gamma, lamlist = lamlist, maxiter = maxiter , eps = eps, penalty = penalty)$path

  # evaluate this on left-out fold:
  x_te <- X[folds[[i]], ]
  x_te <- scale(x_te, center = meanx, scale = FALSE)
  S_te <- crossprod(x_te) / (dim(x_te)[1])

  for (j in seq(nlam)) {
    errs_fit[j, i] <- likelihood(crossprod(path_fit[, , j]), P %*% S_te %*% t(P))
  }
  }

  m_fit <- rowMeans(errs_fit)
  se_fit <- apply(errs_fit, 1, sd) / sqrt(nfolds)
  ibest_fit <- which.min(m_fit)
  i1se_fit <- min(which(m_fit < m_fit[ibest_fit] + se_fit[ibest_fit]))

  fit_cv <- relaxchol(X = X, P = P, lambda = theta[ibest_fit, 1], gamma = theta[ibest_fit, 2],
                  initL = path_fit[, , ibest_fit], maxiter = maxiter, eps = eps, penalty = penalty)


  return(list(errs_fit = errs_fit, folds = folds, lamlist = lamlist, lambda = theta[ibest_fit, 1], gamma = theta[ibest_fit, 2],
            ibest_fit = ibest_fit, i1se_fit = i1se_fit, L_fit = fit_cv))
}

makefolds <- function(n, nfolds) {
  # divides the indices 1:n into nfolds random folds of about the same size.
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}

likelihood <- function (Omega, S){
  # Calculate the negative log-Gaussian likelihood with
  # precision matrix Omega and sample covariance S
  return(-determinant(Omega, logarithm = TRUE)$modulus[1] + sum(S*Omega))
}



################################################
######### Example ##########################
# require(varband)
# set.seed(123)
# n <- 50
# true <- varband_gen(p = 15, block = 5)
# x <- sample_gen(L = true, n = n)
# S <- crossprod(scale(x, center = TRUE, scale = FALSE)) / n
# # init <- diag(1/sqrt(diag(S)))
# # # unweighted estimate
# # L_unweighted <- varband(S, lambda = 0.1, init, w = FALSE)
# # # weighted estimate
# # L_weighted <- varband(S, lambda = 0.1, init, w = TRUE)
# # # lasso estimate
# # L_lasso <- varband(S, lambda = 0.1, init, w = TRUE, lasso = TRUE)
# # ## MCP
# dbug(relaxchol_CV)
# P = diag(1, nrow = 15)
# L_mcp = relaxchol_CV(X = x, gamma = 2, P = diag(1, nrow = 15), lamlist = NULL, nlam = 60,folds = NULL, nfolds = 5, maxiter = 100, penalty = "MCP", eps = 1e-4)

# L_lasso = relaxchol_CV(X = x, gamma = 2, P = diag(1, nrow = 15), lamlist = NULL, nlam = 60,folds = NULL, nfolds = 5, maxiter = 100, penalty = "lasso", eps = 1e-4)
#
# matimage(true)
# matimage(L_lasso$L_fit)
#
# relaxchol_path(X = x, P, gamma = 2, maxiter = 100, eps = 1e-4, penalty = "MCP",
#                          lamlist = NULL, nlam = 60, flmin = 0.01)
