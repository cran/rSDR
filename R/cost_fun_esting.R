cost_fun_esting<-function(C_value,X_train,X_test,Y_test,alpha){
  # Pairwise euclidean distance nxn for response Y
  PDY <- as.matrix(stats::dist(Y_test, method = "euclidean"))
  #update the power of euclidean distance
  KY <- PDY^alpha  # Assuming alpha = 1 for simplicity
  B <- KY - rowMeans(KY) - colMeans(KY) + mean(KY) #B_kl in equation (3)
  N <- stats::cov(X_train)
  N1 <- expm::sqrtm(N)

  # Standardize X
  #Inv.sigma<-solve(N1)
  Z <- X_test %*% solve(N1) #Z is n x d, X is nxp, solve(N1) is pxp.
  ZC <- Z %*% C_value
  PDZC <- as.matrix(stats::dist(ZC))
  PDZC <- PDZC + diag(1, nrow(PDZC))
  f <- -mean((PDZC^alpha) * B)
  f
}
