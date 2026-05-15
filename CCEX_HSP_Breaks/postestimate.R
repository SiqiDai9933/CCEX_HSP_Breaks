library(Matrix)
library(MASS)

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
