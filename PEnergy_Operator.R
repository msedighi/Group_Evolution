PEnergy_Operator <- function(r, f, ...) {
  l <- length(r)
  
  #Defining the difference matrix
  diff_matrix <- t(matrix(rep(r,l),l,l)) - matrix(rep(r,l),l,l)
  
  #Defining the force matrix
  out <- f(Mod(diff_matrix))
  #Make the self force zero
  diag(out) <- 0
  
  out
}