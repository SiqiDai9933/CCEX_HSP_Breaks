library(Matrix)
library(MASS)


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
