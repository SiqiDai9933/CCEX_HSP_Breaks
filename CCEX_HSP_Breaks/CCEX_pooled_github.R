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

run_dp <- function(ssr_mat, m_max, T, min_len = q) {
  opt_ssr <- matrix(Inf, m_max + 1, T); pos_mat <- matrix(0, m_max + 1, T)
  for (n in min_len:T) opt_ssr[1, n] <- ssr_mat[1, n]
  for (m in 2:(m_max + 1)) {
    for (n in (m * min_len):T) {
      j_range <- ((m - 1) * min_len):(n - min_len)
      cand <- sapply(j_range, function(j) opt_ssr[m - 1, j] + ssr_mat[j + 1, n])
      opt_ssr[m, n] <- min(cand); pos_mat[m, n] <- j_range[which.min(cand)] 
    }
  }
  return(list(opt_ssr = opt_ssr, pos_mat = pos_mat))
}


backtrack <- function(pos_mat, m_opt, T) {
  breaks <- c()
  curr_t <- T
  curr_m <- m_opt + 1 
  while (curr_m > 1) {
    j <- pos_mat[curr_m, curr_t]
    breaks <- c(j, breaks)
    curr_t <- j
    curr_m <- curr_m - 1
  }
  return(breaks)
}

postestimate <- function(XT, YT, N, T, minTb) {
  X <- matrix(0, N * T, k); Y <- matrix(0, N * T, 1); X_bar <- matrix(0, T, k); Y_bar <- matrix(0, T, 1)
  for(i in 1:N){
    select <- rep(0, T)
    for(t in 1:T){
      select[t] <- (t - 1) * N + i
      X_bar[t, ] <- colMeans(XT[((t - 1) * N + 1):(t * N), ])
      Y_bar[t, 1] <- mean(YT[((t - 1) * N + 1):(t * N), ])
    }
    X[((i - 1) * T + 1):(i * T), ] <- XT[select, ]
    Y[((i - 1) * T + 1):(i * T), ] <- YT[select, ]
  }
  I_T <- diag(T); M_X_bar <- I_T - tcrossprod(X_bar %*% ginv(crossprod(X_bar, X_bar)), X_bar)
  SW <- genSW(h, N)
  m_est <- length(minTb); Theta0 <- matrix(0, N, (m_est + 1) * (k + 1))
  
  for (i in (1:N)) {
    Qi <- matrix(kronecker(I_T, matrix(SW[i, ], 1, N)) %*% XT, T, k); Qi <- cbind(X[(T*(i-1)+1):(i*T), ], Qi)
    WY_i <- matrix(kronecker(I_T, matrix(SW[i, ], 1, N)) %*% YT, T, 1); Xast_i <- cbind(WY_i, X[(T*(i-1)+1):(i*T), ])
    
    if (m_est == 1) {
      T1 <- minTb; T2 <- T - minTb
      X_bar_1 <- matrix(0, T1, k); X_bar_2 <- matrix(0, T2, k)
      for(t in 1:T1) X_bar_1[t,] <- colMeans(XT[((t-1)*N+1):(t*N),])
      for(t in (T1+1):T) X_bar_2[(t-T1),] <- colMeans(XT[((t-1)*N+1):(t*N),])
      M_X_bar_1 <- diag(T1) - tcrossprod(X_bar_1 %*% ginv(crossprod(X_bar_1, X_bar_1)), X_bar_1)
      M_X_bar_2 <- diag(T2) - tcrossprod(X_bar_2 %*% ginv(crossprod(X_bar_2, X_bar_2)), X_bar_2)
      Qi_1 <- as.matrix(M_X_bar_1 %*% Qi[1:minTb, ]); Qi_2 <- as.matrix(M_X_bar_2 %*% Qi[(minTb + 1):T, ])
      Xast_hat_i_1 <- Qi_1 %*% (ginv(crossprod(Qi_1, Qi_1))) %*% t(Qi_1) %*% Xast_i[1:minTb, ]
      Xast_hat_i_2 <- Qi_2 %*% (ginv(crossprod(Qi_2, Qi_2))) %*% t(Qi_2) %*% Xast_i[(minTb + 1):T, ]
      Theta_i_1 <- ginv(crossprod(Xast_hat_i_1)) %*% crossprod(Xast_hat_i_1, Y[((i-1)*T+1):((i-1)*T+minTb), ])
      Theta_i_2 <- ginv(crossprod(Xast_hat_i_2)) %*% crossprod(Xast_hat_i_2, Y[((i-1)*T+minTb+1):(i*T), ])
      Theta0[i, ] <- cbind(t(Theta_i_1), t(Theta_i_2))
    }
    if (m_est == 2) {
      minTb1 <- minTb[1]; minTb2 <- minTb[2]; T1 <- minTb1; T2 <- minTb2 - minTb1; T3 <- T - minTb2
      X_bar_1 <- matrix(0, T1, k); X_bar_2 <- matrix(0, T2, k); X_bar_3 <- matrix(0, T3, k)
      for(t in 1:T1) X_bar_1[t,] <- colMeans(XT[((t-1)*N+1):(t*N),])
      for(t in (T1+1):(T1+T2)) X_bar_2[(t-T1),] <- colMeans(XT[((t-1)*N+1):(t*N),])
      for(t in (T1+T2+1):T) X_bar_3[(t-T1-T2),] <- colMeans(XT[((t-1)*N+1):(t*N),])
      M_X_bar_1 <- diag(T1) - tcrossprod(X_bar_1 %*% ginv(crossprod(X_bar_1, X_bar_1)), X_bar_1)
      M_X_bar_2 <- diag(T2) - tcrossprod(X_bar_2 %*% ginv(crossprod(X_bar_2, X_bar_2)), X_bar_2)
      M_X_bar_3 <- diag(T3) - tcrossprod(X_bar_3 %*% ginv(crossprod(X_bar_3, X_bar_3)), X_bar_3)
      Qi_1 <- as.matrix(M_X_bar_1 %*% Qi[1:minTb1, ]); Qi_2 <- as.matrix(M_X_bar_2 %*% Qi[(minTb1+1):minTb2, ]); Qi_3 <- as.matrix(M_X_bar_3 %*% Qi[(minTb2+1):T, ])
      Xast_hat_i_1 <- Qi_1 %*% (ginv(crossprod(Qi_1, Qi_1))) %*% t(Qi_1) %*% Xast_i[1:minTb1, ]
      Xast_hat_i_2 <- Qi_2 %*% (ginv(crossprod(Qi_2, Qi_2))) %*% t(Qi_2) %*% Xast_i[(minTb1+1):minTb2, ]
      Xast_hat_i_3 <- Qi_3 %*% (ginv(crossprod(Qi_3, Qi_3))) %*% t(Qi_3) %*% Xast_i[(minTb2+1):T, ]
      Theta_i_1 <- ginv(crossprod(Xast_hat_i_1)) %*% crossprod(Xast_hat_i_1, Y[((i-1)*T+1):((i-1)*T+minTb1), ])
      Theta_i_2 <- ginv(crossprod(Xast_hat_i_2)) %*% crossprod(Xast_hat_i_2, Y[((i-1)*T+minTb1+1):((i-1)*T+minTb2), ])
      Theta_i_3 <- ginv(crossprod(Xast_hat_i_3)) %*% crossprod(Xast_hat_i_3, Y[((i-1)*T+minTb2+1):(i*T), ])
      Theta0[i, ] <- cbind(t(Theta_i_1), t(Theta_i_2), t(Theta_i_3))
    }
  }
  Theta0
}


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


# max m_hat=3 
#input:
# NT \times k dimentional XT=(X'_{.1},...,X'_{.T})';
# NT \times 1 dimentional YT=(Y'_{.1},...,Y'_{.T})';
# N: number of cross sections;
# T: number of time periods;
# k: number of regressors.

#output:
# m_hat: estimate of number of breaks; 
# (we set the maximum number of m_hat is 3 here, it can extend to any numbers following the same logic of our code.)
# minTb_out: estimates of the dates of breaks;
# delta_hat: estimates of parameters in each regime

CCEX_HSP_Breaks<- function(YT,XT,N,T,k) {

  q <- k + 1 
  
  ssr_mat <- get_ssr_matrix(YT, XT, N, T, q=q)
  dp_res <- run_dp(ssr_mat, 3, T) 
  
  best_breaks_list <- list()
  for (m_idx in 1:3) {
    best_breaks_list[[m_idx]] <- backtrack(dp_res$pos_mat, m_idx, T)
  }
  

  f_1_0<- calc_SupF_Pooled(YT,XT, N, T, q, 1, best_breaks_list[[1]])
  f_2_1 <- calc_F_seq_Pooled(YT,XT, N, T, q, 1, best_breaks_list[[1]])
  f_3_2 <- calc_F_seq_Pooled(YT,XT, N, T, q, 2, best_breaks_list[[2]])
  
  
  m_hat <- 0
  
  if (f_1_0 > 13.98) { 
    if (f_2_1 > 15.72) {
      if (f_3_2 > 16.83) m_hat <- 3 else m_hat <- 2
    } else {
      m_hat <- 1
    }
  }
  
  if (m_hat > 0) {
    num_to_fill <- min(m_hat, 2)
    minTb_out[1:num_to_fill] <- best_breaks_list[[m_hat]][1:num_to_fill]
  }
  
    minTb_est <- best_breaks_list[[m_hat]]
    delta_hat <- postestimate(XT, YT, N, T, minTb_est)

  
  results<-cbind(m_hat, minTb_out, delta_hat)
  results
}