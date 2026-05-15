library(Matrix)
library(MASS)


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