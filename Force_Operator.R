Force_Operator <- function(r, f, ...) {
  l <- length(r)
  #Defining the difference matrix
  diff_matrix <- t(matrix(rep(r,l),l,l)) - matrix(rep(r,l),l,l)
  #Defining the force matrix
  f_matrix <- f(Mod(diff_matrix))/Mod(diff_matrix)
  out <- f_matrix*diff_matrix
  #Make the self force zero
  diag(out) <- 0
  
  out
}