library(Matrix)
library(MASS)


get_ssr_matrix <- function(YT, XT, N, T, q) {
  ssr_mat <- matrix(Inf, T, T)
  
  X_bar <- matrix(0, T, k);Y_bar<-matrix(0,T,1)
  for (t in 1:T) {
    X_bar[t, ] <- colMeans(XT[((t - 1) * N + 1):(t * N), ])
    Y_bar[t,1]<-mean(YT[((t-1)*N+1):(t*N),])
  }
  
  
  M_X_bar<-diag(T)-tcrossprod(X_bar%*%ginv(crossprod(X_bar,X_bar)),X_bar)#CCEX
  M_X_bar_NT<-bdiag(replicate(N,M_X_bar,simplify=FALSE))
  
  X <- matrix(0, N*T, k); Y <- matrix(0, N*T, 1); WY <- matrix(0, N*T, 1)
  SW <- genSW(h, N)
  for (i in 1:N) {
    idx_T <- seq(i, (T - 1) * N + i, by = N)
    WY_i <- kronecker(diag(T), matrix(SW[i, ], 1, N)) %*% YT
    X[((i - 1) * T + 1):(i * T), ] <- XT[idx_T, ]#X_1.....,X__N.
    Y[((i - 1) * T + 1):(i * T), ] <-  YT[idx_T, ]
    WY[((i - 1) * T + 1):(i * T), ] <-  WY_i
  }
  X_global<-cbind(X,WY)
  Y_global<-Y
  
  min_obs <- q
  for (i in 1:(T - min_obs + 1)) {
    for (j in (i + min_obs - 1):T) {
      curr_ssr <- 0
      for (u in 1:N) {
        idx <- ((u - 1) * T + i):((u - 1) * T + j)
        X_seg <- matrix(0, N * T, q);Y_seg <- matrix(0, N * T, 1)
        X_seg[idx, ]<- X_global[idx, ]#extend to T*q
        Y_seg[idx, ]<- Y_global[idx, ]
        yu <- M_X_bar%*% Y_seg[((u - 1) * T+1):(u*T), ]
        xu <- M_X_bar%*% X_seg[((u - 1) * T+1):(u*T), ]
        curr_ssr <- curr_ssr + sum((yu - xu %*% ginv(t(xu) %*% xu) %*% t(xu) %*% yu)^2)
      }
      ssr_mat[i, j] <- curr_ssr
    }
  }
  return(ssr_mat)
}