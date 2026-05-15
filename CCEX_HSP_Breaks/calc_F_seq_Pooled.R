library(Matrix)
library(MASS)


calc_F_seq_Pooled <- function(YT, XT, N, T, q, m, current_brks) {
  
  brks_full_m <- c(0, sort(current_brks), T)
  n_regimes_m <- length(brks_full_m) - 1
  
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
  
  
  global_min_ssr <- Inf
  global_best_tau <- NA
  
  for (j in 1:n_regimes_m) {
    t_start <- brks_full_m[j] + 1
    t_end   <- brks_full_m[j + 1]
    Tj_len  <- t_end - t_start + 1
    
    h_trim <- max(floor(0.15 * Tj_len), q + 2)
    if (Tj_len < (2 * h_trim)) next
    
    search_range <- (t_start + h_trim - 1):(t_end - h_trim)
    
    for (tau in search_range) {
      temp_brks <- sort(c(current_brks, tau))
      temp_brks_full <- c(0, temp_brks, T)
      n_regimes_new <- m + 2
      total_params_new <- n_regimes_new * q
      
      X_seg <- matrix(0, nrow = N * T, ncol = total_params_new)
      for (r in 1:n_regimes_new) {
        tr_start <- temp_brks_full[r] + 1
        tr_end   <- temp_brks_full[r + 1]
        col_idx <- ((r - 1) * q + 1):(r * q)
        for (i in 1:N) {
          rows_i_r <- (i - 1) * T + (tr_start:tr_end)
          X_seg[rows_i_r, col_idx] <- X_global[rows_i_r, ]
        }
      }
      
      X_tilde <- as.matrix(M_X_bar_NT %*% X_seg)
      Y_tilde <- M_X_bar_NT %*% Y_global
      
      XtX_inv <- ginv(t(X_tilde) %*% X_tilde)
      P_X <- X_tilde %*% XtX_inv %*% t(X_tilde)
      resid <- Y_tilde - P_X %*% Y_tilde
      ssr <- sum(resid^2)
      
      if (ssr < global_min_ssr) {
        global_min_ssr <- ssr
        global_best_tau <- tau
      }
    }
  }
  
  if (is.na(global_best_tau)) {
    return(NA)
  }
  
  
  temp_brks <- sort(c(current_brks, global_best_tau))
  temp_brks_full <- c(0, temp_brks, T)
  n_regimes_new <- m + 2
  total_params <- n_regimes_new * q
  
  X_seg <- matrix(0, nrow = N * T, ncol = total_params)
  for (r in 1:n_regimes_new) {
    tr_start <- temp_brks_full[r] + 1
    tr_end   <- temp_brks_full[r + 1]
    col_idx <- ((r - 1) * q + 1):(r * q)
    for (i in 1:N) {
      rows_i_r <- (i - 1) * T + (tr_start:tr_end)
      X_seg[rows_i_r, col_idx] <- X_global[rows_i_r, ]
    }
  }
  
  X_tilde <- as.matrix(M_X_bar_NT %*% X_seg)
  Y_tilde <- M_X_bar_NT %*% Y_global
  
  
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
    
    score_sum <- score_sum + score_i %*% t(score_i)
  }
  
  Phi_hat <- score_sum / (N * T)
  
  Omega_hat <- XtX / (N * T)
  
  Omega_inv <- ginv(Omega_hat)
  V_theta_pooled <- Omega_inv %*% Phi_hat %*% Omega_inv
  
  
  new_break_pos <- which(temp_brks == global_best_tau)
  R_mat <- matrix(0, q, n_regimes_new * q)
  R_mat[, ((new_break_pos - 1) * q + 1):(new_break_pos * q)] <- diag(q)
  R_mat[, (new_break_pos * q + 1):((new_break_pos + 1) * q)] <- -diag(q)
  
  delta_diff <- R_mat %*% theta_pooled_stack
  V_delta <- R_mat %*% V_theta_pooled %*% t(R_mat)
  
  wald_val <- as.numeric(t(delta_diff) %*% ginv(V_delta) %*% delta_diff)
  
  prefix <- (N * (T - (m + 2) * q) - (m + 2) * q) / q
  F_stat <- prefix * wald_val
  return(F_stat)
}
