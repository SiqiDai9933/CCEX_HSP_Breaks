library(Matrix)
library(MASS)

calc_SupF_Pooled <- function(YT, XT, N, T, q, m, brks) {
  
  brks_full <- c(0, sort(brks), T)
  n_regimes <- m + 1
  total_params <- n_regimes * q
  
  X_bar <- matrix(0, T, k)
  Y_bar <- matrix(0, T, 1)
  for (t in 1:T) {
    X_bar[t, ] <- colMeans(XT[((t - 1) * N + 1):(t * N), ])
    Y_bar[t, 1] <- mean(YT[((t - 1) * N + 1):(t * N), ])
  }
  
  M_X_bar <- diag(T) - tcrossprod(X_bar %*% ginv(crossprod(X_bar, X_bar)), X_bar)
  M_X_bar_NT <- bdiag(replicate(N, M_X_bar, simplify = FALSE))
  
  
  X <- matrix(0, N * T, k)
  Y <- matrix(0, N * T, 1)
  WY <- matrix(0, N * T, 1)
  SW <- genSW(h, N)
  
  for (i in 1:N) {
    idx_T <- seq(i, (T - 1) * N + i, by = N)
    WY_i <- kronecker(diag(T), matrix(SW[i, ], 1, N)) %*% YT
    X[((i - 1) * T + 1):(i * T), ] <- XT[idx_T, ]
    Y[((i - 1) * T + 1):(i * T), ] <- YT[idx_T, ]
    WY[((i - 1) * T + 1):(i * T), ] <- WY_i
  }
  X_global <- cbind(WY, X)
  Y_global <- Y
  
  Y_tilde <- M_X_bar_NT %*% Y_global
  
  X_seg <- matrix(0, nrow = N * T, ncol = total_params)
  for (j in 1:n_regimes) {
    t_start <- brks_full[j] + 1
    t_end   <- brks_full[j + 1]
    col_idx <- ((j - 1) * q + 1):(j * q)
    for (i in 1:N) {
      rows_i_j <- (i - 1) * T + (t_start:t_end)
      X_seg[rows_i_j, col_idx] <- X_global[rows_i_j, ]
    }
  }
  
  X_tilde <- as.matrix(M_X_bar_NT %*% X_seg)
  
  
  XtX <- t(X_tilde) %*% X_tilde
  XtY <- t(X_tilde) %*% Y_tilde
  theta_pooled_stack <- as.numeric(ginv(XtX) %*% XtY)
  
  residuals <- as.numeric(Y_tilde - X_tilde %*% theta_pooled_stack)
  
  score_sum <- matrix(0, total_params, total_params)
  
  for (i in 1:N) {
    
    idx_range <- ((i - 1) * T + 1):(i * T)
    
    
    X_i <- X_tilde[idx_range, ]  # T × total_params
    e_i <- residuals[idx_range]   # T × 1
    
    
    score_i <- t(X_i) %*% e_i
    
    # sum up score_i * score_i'
    score_sum <- score_sum + score_i %*% t(score_i)
  }
  
  # Phi_hat
  Phi_hat <- score_sum / (N * T)
  
  Omega_hat <- XtX / (N * T)
  
  Omega_inv <- ginv(Omega_hat)
  V_theta_pooled <- Omega_inv %*% Phi_hat %*% Omega_inv
  
  R <- matrix(0, m * q, total_params)
  for (j in 1:m) {
    row_start <- (j - 1) * q + 1
    row_end   <- j * q
    
    col1_start <- (j - 1) * q + 1
    col1_end   <- j * q
    col2_start <- j * q + 1
    col2_end   <- (j + 1) * q
    
    R[row_start:row_end, col1_start:col1_end] <- diag(q)
    R[row_start:row_end, col2_start:col2_end] <- -diag(q)
  }
  
  delta_pooled <- R %*% theta_pooled_stack
  V_delta <- R %*% V_theta_pooled %*% t(R)
  
  wald <- as.numeric(t(delta_pooled) %*% ginv(V_delta) %*% delta_pooled)
  
  df_adj <- (N * (T - (m + 1) * q) - (m + 1) * q) / (m * q)
  F_stat <- df_adj * wald
  return(F_stat)
}
