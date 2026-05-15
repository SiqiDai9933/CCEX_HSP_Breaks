#input:
# NT \times k dimentional XT=(X'_{.1},...,X'_{.T})';
# NT \times 1 dimentional YT=(Y'_{.1},...,Y'_{.T})';
# N: number of cross sections;
# T: number of time periods;
# k: number of regressors.

#output:
# m_hat: estimate of number of breaks; 
# (we set the maximum number of m_hat (m_max) is 3 here, it can extend to any numbers following the same logic of our code.)
# minTb_out: estimates of the dates of breaks;
# delta_hat: estimates of parameters in each regime

library(Matrix)
library(MASS)

CCEX_HSP_Breaks<- function(YT,XT,N,T,k) {
  source("get_ssr_matrix.R")
  source("run_dp.R") 
  source("backtrack.R") 
  source("postestimate.R")
  source("cal_SupF_Pooled.R") 
  source("cal_F_seq_Pooled.R") 
  
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