#' Robust Sufficient Dimension Reduction
#' @description
#' Robust Sufficient Dimension Reduction with alpha-Distance Covariance and Stiefel Manifold Learning for supervised dimension reduction.
#' @usage rSDR(X, Y, d, alpha=0.5,maxiter=1000,tol=1e-7)
#' @param X an \eqn{n \times p} numeric matrix, where \eqn{n} is the number of observations and \eqn{p} is the number of variable.
#' @param Y an \eqn{n \times k} numeric response matrix, where \eqn{k (\geq 1)} is the number of response variables.
#' @param d the number of reduced dimension.
#' @param alpha this parameter represents the exponent applied to the Euclidean distance in the computation of distance covariance. When alpha=1, it corresponds to the classical distance covariance. When 0 < alpha < 1, it is a more robust version by reducing the influence of large values in the distance matrices.
#' @param maxiter maxiter is the maximum number of iterations allowed for the solver (a non-negative integer). See the Max_Iteration parameter in \code{\link[ManifoldOptim]{get.solver.params}} for details.
#' @param tol tol is used to assess convergence, see the Tolerance parameter in \code{\link[ManifoldOptim]{get.solver.params}} for details.
#'
#'
#' @return The returned value is an object of class "rSDR", containing the following components:
#'  \describe{
#'   \item{projected_data}{an \eqn{n \times d} matrix representing the projected data using the rSDR method.}
#'   \item{beta}{a \eqn{p \times d} matrix.Solve \eqn{\boldsymbol{\beta}} by \eqn{\textbf{C}=\Sigma_x^{1/2} \boldsymbol{\beta}}}
#'   \item{C_value}{an optimal of C is obtained by maximizing the target function using ManifoldOptim method.}
#'   \item{f_value}{The value of cost function f is defined as the negative of the target function.}
#' }
#'
#' @export
#' @references
#'  Hsin-Hsiung Huang, Feng Yu & Teng Zhang (19 Feb 2024): Robust sufficient dimension reduction via alpha-distance covariance, Journal of Nonparametric Statistics, DOI:10.1080/10485252.2024.2313137
#'
#' @importFrom stats cov dist
#' @examples
#'
#' library(ManifoldOptim)
#' library(rSDR)
#' utils::data("ionosphere", package = "fdm2id")
#' X<-as.matrix(ionosphere[,c(1:33)])
#' Y<-ifelse(ionosphere[,34]=='b',0,1)
#' Y<-matrix(Y,length(Y),1)
#' set.seed(2435)
#' \donttest{
#' sdr_result<-rSDR(X=X, Y=Y, d=3, alpha=0.3,maxiter=1000,tol=1e-7)
#' }
rSDR <- function(X, Y, d, alpha=0.5,maxiter=1000,tol=1e-7) {

  if (!is.matrix(X)){
    stop('X should be a matrix.')
  }

  if (!(is.matrix(Y)&is.numeric(Y))){
    stop('Y should be a matrix.')
  }

  # Extract dimensions
  n <- nrow(X)
  p <- ncol(X)

  # Initial matrix C for the problem
  # init is of size p x d
  init <- rstiefel::rustiefel(p, d)  # Generate a random orthonormal matrix (C'C=CC'=I) for initialization

  # Ensure n >= p
  if (n < p) {
    stop("The dimension n must be larger than the dimension p.")
  }

  # Pairwise euclidean distance nxn for response Y
  PDY <- as.matrix(stats::dist(Y, method = "euclidean"))
  #update the power of euclidean distance
  KY <- PDY^alpha  # alpha = 1 for simplicity
  B <- KY - rowMeans(KY) - colMeans(KY) + mean(KY) #B_kl in equation (3)
  N <- stats::cov(X)
  N1 <- expm::sqrtm(N)

  # Standardize X
  #Inv.sigma<-solve(N1)
  Z <- X %*% solve(N1) #Z is n x d, X is nxp, solve(N1) is pxp.

  # Define cost function and gradient
  # Cost function
  cost_store <- function(C, Z, B, alpha, store = list()) {
    ZC <- Z %*% C
    #ZC pairwise distance
    PDZC <- as.matrix(stats::dist(ZC))
    PDZC <- PDZC + diag(1, nrow(PDZC))
    f <- -mean((PDZC^alpha) * B)
    store$PDZC <- PDZC
    return(list(f = f, store = store))
  }

  #gradient
  eGrad_store <- function(C, Z, B, alpha, store) {
    PDZC <- store$PDZC
    coef <- alpha * (PDZC^2 + 1e-10)^((alpha - 2) / 2) * B
    G <- -2 * t(Z) %*% ((diag(colSums(coef)) - coef) %*% Z %*% C) / n^2
    return(list(G = G, store = store))
  }

  # Modify the cost and gradient functions to accept vector inputs
  f <- function(c_vec) {
    C <- matrix(c_vec, nrow = p, ncol = d)
    cost_store(C, Z, B, alpha)$f  # Assuming alpha = 1
  }

  g <- function(c_vec) {
    C <- matrix(c_vec, nrow = p, ncol = d)
    G <- eGrad_store(C, Z, B, alpha, cost_store(C, Z, B, alpha)$store)$G
    as.numeric(G)  # Flatten the gradient matrix to a vector
  }

  # Create the problem object using RProblem class
  mod <- Rcpp::Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob <- methods::new(mod$RProblem, f, g)

  # Initial point for optimization (flattened to a vector)
  x0 <- as.numeric(init)  # initial point

  # Define the manifold using stiefel_factory (if necessary)
  mani.defn <- ManifoldOptim::get.stiefel.defn(p, d)

  # Set solver parameters
  solver.params <- ManifoldOptim::get.solver.params()
  solver.params$Max_Iteration <- maxiter
  solver.params$Tolerance <- tol

  # Perform optimization
  result <- ManifoldOptim::manifold.optim(
    prob = prob,
    mani.defn = mani.defn,
    method = "RCG", # RCG for Riemannian Conjugate Gradient
    x0 = x0,
    solver.params = solver.params
  )

  # Extract C (pxd) from results
  C_value <- matrix(result$xopt, nrow = p, ncol = d)
  #Extract a value of f from results
  #f_value<-f(C_value)
  f_value<-result$fval

  ##Since C=sigma_x^(0.5)*beta, beta=sigma_x^(-0.5)%*%C from equation (4)
  beta <- solve(N1) %*% C_value

  projected_data <- X %*% beta  # This projects X onto the estimated SDR directions
  projected_data<-as.data.frame(projected_data)

  # Return the result
  return(list(projected_data=projected_data,beta = beta,C_value=C_value,f_value=f_value))

}
