

#' @title Rank selection via eigenvalue-gap methods
#'
#' @description
#' Select the number of significant principal components of a GMF model via
#' exploitation of eigenvalue-gap methods
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row-specific fixed effects (\eqn{n \times p})
#' @param Z matrix of column-specific fixed effects (\eqn{q \times m})
#' @param maxcomp maximum number of eigenvalues to compute
#' @param family a family as in the \code{\link{glm}} interface (default \code{gaussian()})
#' @param weights matrix of optional weights (\eqn{n \times m})
#' @param offset matrix of optional offsets (\eqn{n \times m})
#' @param method rank selection method
#' @param type.reg regression method to be used to profile out the covariate effects
#' @param type.res residual type to be decomposed
#' @param normalize if \code{TRUE}, standardize column-by-column the residual matrix
#' @param maxiter maximum number of iterations
#' @param parallel if \code{TRUE}, allows for parallel computing using \code{foreach}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#' @param return.eta if \code{TRUE}, return the linear predictor martix
#' @param return.mu if \code{TRUE}, return the fitted value martix
#' @param return.res if \code{TRUE}, return the residual matrix
#' @param return.cov if \code{TRUE}, return the covariance matrix of the residuals
#'
#' @return
#' A list containing the \code{method}, the selected latent rank \code{ncomp},
#' and the eigenvalues used to select the latent rank \code{lambdas}.
#' Additionally, if required, in the output list will also provide the linear predictor
#' \code{eta}, the predicted mean matrix \code{mu}, the residual matrix \code{res}, and
#' the implied residual covariance matrix \code{covmat}.
#'
#' @references
#' Onatski, A. (2010).
#' \emph{Determining the number of factors from empirical distribution of eigenvalues.}
#' Review of Economics and Statistics, 92(4): 1004-1016
#'
#' Ahn, S.C., Horenstein, A.R. (2013).
#' \emph{Eigenvalue ratio test for the number of factors.}
#' Econometrica, 81, 1203-1227
#'
#' Gavish, M., Donoho, D.L. (2014)
#' \emph{The optimal hard thresholding for singular values is 4/sqrt(3).}
#' IEEE Transactions on Information Theory, 60(8): 5040--5053
#'
#' Fan, J., Guo, J. and Zheng, S. (2020).
#' \emph{Estimating number of factors by adjusted eigenvalues thresholding.}
#' Journal of the American Statistical Association, 117(538): 852--861
#'
#' Wang, L. and Carvalho, L. (2023).
#' \emph{Deviance matrix factorization.}
#' Electronic Journal of Statistics, 17(2): 3762-3810
#'
#' @examples
#' library(sgdGMF)
#'
#' # Set the data dimensions
#' n = 100; m = 20; d = 5
#'
#' # Generate data using Poisson, Binomial and Gamma models
#' data_pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
#' data_bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
#' data_gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"), dispersion = 0.25)
#'
#' # Initialize the GMF parameters assuming 3 latent factors
#' ncomp_pois = sgdgmf.rank(data_pois$Y, family = poisson(), normalize = TRUE)
#' ncomp_bin = sgdgmf.rank(data_bin$Y, family = binomial(), normalize = TRUE)
#' ncomp_gam = sgdgmf.rank(data_gam$Y, family = Gamma(link = "log"), normalize = TRUE)
#'
#' # Get the selected number of components
#' print(paste("Poisson:", ncomp_pois$ncomp))
#' print(paste("Binomial:", ncomp_bin$ncomp))
#' print(paste("Gamma:", ncomp_gam$ncomp))
#'
#' # Plot the screeplot used for the component determination
#' oldpar = par(no.readonly = TRUE)
#' par(mfrow = c(3,1))
#' barplot(ncomp_pois$lambdas, main = "Poisson screeplot")
#' barplot(ncomp_bin$lambdas, main = "Binomial screeplot")
#' barplot(ncomp_gam$lambdas, main = "Gamma screeplot")
#' par(oldpar)
#'
#' @export sgdgmf.rank
sgdgmf.rank = function (
    Y,
    X = NULL,
    Z = NULL,
    maxcomp = ncol(Y),
    family = gaussian(),
    weights = NULL,
    offset = NULL,
    method = c("evr", "onatski", "act", "oht"),
    type.reg = c("ols", "glm"),
    type.res = c("deviance", "pearson", "working", "link"),
    normalize = FALSE,
    maxiter = 10,
    parallel = FALSE,
    nthreads = 1,
    return.eta = FALSE,
    return.mu  = FALSE,
    return.res = FALSE,
    return.cov = FALSE
) {
  # Set the selection method
  method = match.arg(method)
  type.reg = match.arg(type.reg)
  type.res = match.arg(type.res)

  # Set the family-specific data transformation
  family = set.family(family)

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)

  # Fill the missing values
  Y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (x)
  })

  # Compute the transformed data
  gY = family$transform(Y)

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = n, ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = m, ncol = 1)
  if (is.null(weights)) weights = matrix(1, nrow = n, ncol = m)
  if (is.null(offset)) offset = matrix(0, nrow = n, ncol = m)

  # Register and open the connection to the clusters
  clust = NULL
  if (parallel) {
    ncores = parallel::detectCores() - 1
    ncores = max(1, min(nthreads, ncores))
    clust = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(clust)
  }

  # Set the required GLM matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu  = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # Initialize the regression coefficients
  eta[] = offset
  B = switch(type.reg,
    "ols" = ols.fit.coef(gY, X, offset = eta),
    "glm" = vglm.fit.coef(Y, X, family = family, weights = weights,
                          offset = eta, parallel = parallel,
                          nthreads = nthreads, clust = clust))

  eta[] = eta + tcrossprod(X, B)
  A = switch(type.reg,
    "ols" = ols.fit.coef(t(gY), Z, offset = t(eta)),
    "glm" = vglm.fit.coef(t(Y), Z, family = family, weights = t(weights),
                          offset = t(eta), parallel = parallel,
                          nthreads = nthreads, clust = clust))

  # Close the connection to the clusters
  if (parallel) parallel::stopCluster(clust)

  # Initialize the linear predictor and the conditional mean matrix
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)
  res[] = switch(type.res,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (Y - mu) * family$mu.eta(eta) / abs(family$variance(mu)),
    "link" = (gY - eta))

  # Compute the covariance/correlation matrix of the residuals
  covmat = if (normalize) cor(res) else cov(res)

  # Select the optimal rank
  eigengap = switch(method,
    "onatski" = eigengap.onatski(covmat, min(maxcomp, n-6, m-6), maxiter),
    "act" = eigengap.act(covmat, n, maxcomp),
    "oht" = eigengap.oht(covmat, n, maxcomp),
    "evr" = eigengap.evr(covmat, maxcomp))

  # Build the output
  out = list()

  out$method = method
  out$ncomp = eigengap$ncomp
  out$lambdas = eigengap$lambdas

  if (return.eta) out$eta = eta
  if (return.mu) out$mu = mu
  if (return.res) out$residuals = res
  if (return.cov) out$covariance = covmat

  # Return the selected rank
  return (out)
}

#' @title Rank selection via eigenvalue ratio maximization
#'
#' @description
#' Select the number of significant principal components of a matrix via the
#' eigenvalue ratio (EVR) maximization method
#'
#' @param covmat matrix to be decomposed
#' @param maxcomp maximum number of eigenvalues to compute
#'
#' @references
#' Ahn, S.C., Horenstein, A.R. (2013).
#' \emph{Eigenvalue ratio test for the number of factors.}
#' Econometrica, 81, 1203-1227
#'
#' @keywords internal
eigengap.evr = function(covmat, maxcomp = 50, thr = 0.95) {

  # Set the matrix dimension
  m = ncol(covmat)

  # Safety check for the number of maximum components
  maxcomp = floor(maxcomp)
  if (maxcomp < 1 | maxcomp > m) {
    maxcomp = max(1, min(m, maxcomp))
    warning("Rank selection: 'maxcomp' set to default value.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }

  # Compute the spectrum of the covariance matrix of Y
  if (maxcomp < m) {
    lambdas = RSpectra::eigs_sym(covmat, maxcomp)$values
  } else {
    lambdas = eigen(covmat)$values
  }

  # Compute the maximum of the eigenvalue ratio
  # (I exclude some eigenvalues if their do not explain enough variance)
  idx = head(cumsum(lambdas) < thr * sum(lambdas), -1)
  ratio = head(lambdas, -1) / tail(lambdas, -1)
  ncomp = which.max(ratio[idx]) + 1

  #Return the selected rank
  list(ncomp = ncomp, lambdas = lambdas, ratio = ratio)
}

#' @title Rank selection via the Onatski method
#'
#' @description
#' Select the number of significant principal components of a matrix via the
#' Onatski method
#'
#' @param covmat matrix to be decomposed
#' @param maxcomp maximum number of eigenvalues to compute
#' @param maxiter maximum number of iterations
#'
#' @references
#' Onatski, A. (2010).
#' \emph{Determining the number of factors from empirical distribution of eigenvalues.}
#' Review of Economics and Statistics, 92(4): 1004-1016
#'
#' @keywords internal
eigengap.onatski = function (covmat, maxcomp = 50, maxiter = 100) {

  # Set the matrix dimension
  m = ncol(covmat)

  # Safety check for the number of maximum components
  if (maxcomp > m - 5) {
    maxcomp = max(1, min(m - 5, floor(m / 2)))
    warning("Rank selection: 'maxcomp' set to default value.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }

  # Compute the spectrum of the covariance matrix of Y
  # lambdas = eigen(covmat)$values
  if (maxcomp < m - 5) {
    lambdas = RSpectra::eigs_sym(covmat, maxcomp + 5)$values
  } else {
    lambdas = eigen(covmat)$values
  }

  # Initialize the loop parameters
  tol = 1e+03
  iter = 0
  ju = maxcomp + 1
  ncomp = maxcomp

  # Onatski selection loop
  while (tol >= 1) {
    j = ju
    yreg = lambdas[j:(j+4)]
    xreg = (j + (-1:3))^(2/3)
    delta = 2 * abs(stats::cov(yreg, xreg)) / stats::var(xreg)
    flag = which(-diff(lambdas) >= delta)
    ncomp = ifelse(length(flag) == 0, 0, utils::tail(flag[flag <= maxcomp], 1))
    ju = ncomp + 1
    tol = abs(ju - j)
    iter = iter + 1
    if (iter > maxiter) break
  }

  # Check if the search was successful
  success = iter < maxiter
  ncomp = max(1, ncomp)

  # Return the selected rank
  list(ncomp = ncomp, lambdas = lambdas, delta = delta,
       niter = iter, success = success)
}


#' @title Rank selection via adjust correlation thresholding
#'
#' @description
#' Select the number of significant principal components of a matrix via adjust
#' correlation threshold (ACT)
#'
#' @param covmat matrix to be decomposed
#' @param nobs number of observations used to compute the covariance matrix
#' @param maxcomp maximum number of eigenvalues to compute
#'
#' @references
#' Fan, J., Guo, j. and Zheng, S. (2020).
#' \emph{Estimating number of factors by adjusted eigenvalues thresholding.}
#' Journal of the American Statistical Association, 117(538): 852--861
#'
#' @keywords internal
eigengap.act = function (covmat, nobs, maxcomp = NULL) {
  # Set the data dimensions
  n = nobs
  p = ncol(covmat)
  d = ifelse(is.null(maxcomp), p, maxcomp)

  # Convert the covariance matrix to a correlation matrix
  cormat = stats::cov2cor(covmat)

  # Compute the spectrum of the correlation matrix of Y
  if (p == d) {
    lambdas = eigen(cormat)$values
  } else {
    lambdas = RSpectra::eigs_sym(cormat, d)$values
    lambda0 = p - sum(lambdas)
    lambdas = c(lambdas, rep(lambda0, p - d))
  }

  # ACT selection loop
  m = rep(0, times = d-1)
  for (j in 1:(d-1)) {
    delta1 = 1 / (lambdas[(j+1):d] - lambdas[j])
    delta2 = 1 / ((3 * lambdas[j] + lambdas[j+1]) / 4 - lambdas[j])
    m[j] = 1 / (d-j) * (sum(delta1) + delta2)
  }

  # Rank selection
  rho = (d - 1:(d-1)) / (n-1)
  m1 = rho * m - (1 - rho) / lambdas[1:(d-1)]
  adj.lambdas = - 1 / m1
  thr = 1 + sqrt(d / n)
  ncomp = sum(adj.lambdas > thr)
  ncomp = max(1, ncomp)

  # Return the selected rank
  list(ncomp = ncomp, lambdas = lambdas,
       adj.lambdas = adj.lambdas, threshold = thr)
}


#' @title Rank selection via optimal hard thresholding
#'
#' @description
#' Select the number of significant principal components of a matrix via optimal
#' hard thresholding (OHT)
#'
#' @param covmat matrix to be decomposed
#' @param nobs number of observations used to compute the covariance matrix
#' @param maxcomp maximum number of eigenvalues to compute
#'
#' @references
#' Gavish, M., Donoho, D.L. (2014)
#' \emph{The optimal hard thresholding for singular values is 4/sqrt(3).}
#' IEEE Transactions on Information Theory, 60(8): 5040--5053
#'
#' @keywords internal
eigengap.oht = function (covmat, nobs, maxcomp = NULL) {

  # lambdas = eigen(covmat)$values
  # beta = nobs / ncol(covmat)
  # beta = ifelse(beta > 1, beta, 1 / beta)
  # omega = 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43
  # thr = omega * sqrt(median(lambdas))
  # ncomp = max(1, sum(sqrt(lambdas) > thr))

  lambdas = eigen(covmat)$values
  thr = 2.858 * stats::median(lambdas)
  ncomp = max(1, sum(lambdas > thr))

  # Return the selected rank
  list(ncomp = ncomp, lambdas = lambdas, threshold = thr)
}


