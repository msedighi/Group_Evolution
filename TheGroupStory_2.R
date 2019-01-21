source("C:/Users/tpl97/Documents/R/Force_Operator.R")
source("C:/Users/tpl97/Documents/R/PEnergy_Operator.R")

Num_Points <- 4;
Time_range <- 500;
delta_Time <- .01

Position <- array(NA, c(Num_Points,Time_range))
Velocity <- array(NA, c(Num_Points,Time_range))
Force_vector <- array(NA, c(Num_Points,Time_range))

Position_Spectrum <- array(NA, c(Num_Points,Time_range))
Velocity_Spectrum <- array(NA, c(Num_Points,Time_range))
Force_Spectrum <- array(NA, c(Num_Points,Time_range))

Hamiltonian_Change_Spectrum <- array(NA, c(Num_Points,Time_range))
Principal_Fields <- array(NA, c(Num_Points,Num_Points,Time_range))
Right_Principal_Fields <- array(NA, c(Num_Points,Num_Points,Time_range))
Left_Principal_Fields <- array(NA, c(Num_Points,Num_Points,Time_range))
Principal_Rotation <- array(NA, c(Num_Points,Num_Points,Time_range))

#Initial Conditions
t=0
#Position
x <- runif(Num_Points, -50, 50)
y <- runif(Num_Points, -50, 50)
# x <- c(1,1,-1,-1)
# y <- c(1,-1,1,-1)
z <- complex(real = x, imaginary = y)
z <- z - mean(z)
#Velocity
Vx <- runif(Num_Points, -10, 10)
Vy <- runif(Num_Points, -10, 10)
Vz <- complex(real = Vx, imaginary = Vy)
Vz <- Vz - mean(Vz)
#Defining the force function between points
# f = function(r) -r #Spring
# e = function(r) .5*r^2 #Spring

f = function(r) -r^(-2)*100 #Gravity
e = function(r) -r*100 #Gravity

# f = function(r) -1/r*100 #1/r

while(t < Time_range) {

  t = t + 1
  
  PE <- PEnergy_Operator(z,e)  
  F <- Force_Operator(z,f)
  Vac <- matrix(rep(1,Num_Points),Num_Points,1)
  Force_vector[,t] <- -F %*% Vac
  
  #Building "Force Hamiltonian" Operator
  Force_Hamiltonian <- F
  diag(Force_Hamiltonian) <- Force_vector[,t]
  #Finding "Force Hamiltonian" Spectrum
  Force_svd <- svd(Force_Hamiltonian)
  # Force_eigenstructure <- eigen(Force_Hamiltonian, symmetric=FALSE, only.values = FALSE, EISPACK = FALSE)
  
  
  #Building "Energy Change Operator" (Hamiltonian Change)
  Hamiltonian_Change <- diag(Re(Vz)) %*% Re(Force_Hamiltonian) + diag(Im(Vz)) %*% Im(Force_Hamiltonian)
  PEnergy_Change <- Hamiltonian_Change-t(Hamiltonian_Change)
  #Finding "Hamiltonian Change" Spectrum
  Hamiltonian_Change_svd <- svd(Hamiltonian_Change)
  Hamiltonian_Change_eigenstructure <- eigen(Hamiltonian_Change, symmetric=FALSE, only.values = FALSE)
  
  Position[,t] <- z
  Velocity[,t] <- Vz

  # Principal_Fields[,,t] <- Force_svd$v
  # Position_Spectrum[,t] <- diag(Conj(t(Force_svd$u)) %*% diag(Position[,t]) %*% Force_svd$v)
  # Velocity_Spectrum[,t] <- diag(Conj(t(Force_svd$u)) %*% diag(Velocity[,t]) %*% Force_svd$v)
  Force_Spectrum[,t] <- Force_svd$d # These are all real values!
  
  Principal_Fields[,,t] <- Hamiltonian_Change_svd$v
  Left_Principal_Fields[,,t] <- Hamiltonian_Change_svd$u
  Principal_Rotation[,,t] <- t(Hamiltonian_Change_svd$u) %*% Hamiltonian_Change_svd$v
  Position_Spectrum[,t] <- diag(t(Hamiltonian_Change_svd$v) %*% diag(Position[,t]) %*% Hamiltonian_Change_svd$v)
  Velocity_Spectrum[,t] <- diag(t(Hamiltonian_Change_svd$v) %*% diag(Velocity[,t]) %*% Hamiltonian_Change_svd$v)
  
  # Force_Spectrum[,t] <- diag(Conj(t(Hamiltonian_Change_svd$u)) %*% Force_Hamiltonian %*% Hamiltonian_Change_svd$v)
  Hamiltonian_Change_Spectrum[,t] <- Hamiltonian_Change_svd$d # These are all real values

  # Principal_Fields[,,t] <- Force_eigenstructure$vectors
  # Position_Spectrum[,t] <- diag(Conj(t(Force_eigenstructure$vectors)) %*% diag(Position[,t]) %*% Force_eigenstructure$vectors)
  # Velocity_Spectrum[,t] <- diag(Conj(t(Force_eigenstructure$vectors)) %*% diag(Velocity[,t]) %*% Force_eigenstructure$vectors)
  # Force_Spectrum[,t] <- Force_eigenstructure$values

  # Principal_Fields[,,t] <- Hamiltonian_Change_eigenstructure$vectors
  # Position_Spectrum[,t] <- diag(t(Hamiltonian_Change_eigenstructure$vectors) %*% diag(Position[,t]) %*% Hamiltonian_Change_eigenstructure$vectors)
  # Velocity_Spectrum[,t] <- diag(t(Hamiltonian_Change_eigenstructure$vectors) %*% diag(Velocity[,t]) %*% Hamiltonian_Change_eigenstructure$vectors)

  # Hamiltonian_Change_Spectrum[,t] <- Hamiltonian_Change_eigenstructure$values
  
  #Time Evoloution
  z = z + Vz * delta_Time
  Vz = Vz + Force_vector[,t] * delta_Time
  
}

# win.graph()
# plot(Position, col=rainbow(Num_Points))
# title(main="Position")
# win.graph()
# plot(Velocity, col=rainbow(Num_Points))
# title(main="Velocity")
# win.graph()
# plot(Force_vector, col=rainbow(Num_Points))
# title(main="Force")

# win.graph()
# plot(Position_Spectrum)
# title(main="Position Spectrum")
# win.graph()
# plot(Velocity_Spectrum)
# title(main="Velocity Spectrum")
# win.graph()
# For Eigen Structure:
# plot(Force_Spectrum)
# For SVD:
plot(t(matrix(rep(1:Time_range, Num_Points),Time_range,Num_Points)), Force_Spectrum, col=matrix(rep(rainbow(Num_Points),Time_range),Num_Points,Time_range))
title(main="Force Spectrum")

if (TRUE) {
  win.graph()
  plot(Position, col=rainbow(Num_Points))
  points(Position_Spectrum)
  title(main="Position")
  win.graph()
  plot(Velocity, col=rainbow(Num_Points))
  points(Velocity_Spectrum)
  title(main="Velocity")
  win.graph()
  plot(Force_vector, col=rainbow(Num_Points))
  points(Force_Spectrum)
  title(main="Force")
}

# Problems I am interested in:
#   1-Studying eigenvectors
#   2-Differentiating the eigenvalies based on eigenvectors. Requires following eigenvectors in time (?)
#   3-Plotting observables in time (third axis) or animation
#   4-Defining an "Energy Operator" & finding its relationship with other operators specially the "Force Operator"
#   5-Studying the Lee Algebra of operators, specially the commutators
#   6-Studying the Wigner diagrams (if possible for 1-D particles to be visualized easily)
#   7-The problem of scale and all its interesting following results
#   8-The generalization of the problem to 3D (quaternions? exterior vector product?)

win.graph()
plot(t(matrix(rep(1:Time_range, Num_Points),Time_range,Num_Points)), Hamiltonian_Change_Spectrum, col=matrix(rep(rainbow(Num_Points),Time_range),Num_Points,Time_range))
title(main="Hamiltonian Change Spectrum")

win.graph()
plot(t(matrix(rep(1:Time_range, Num_Points),Time_range,Num_Points)), Principal_Fields[,2,], col=matrix(rep(rainbow(Num_Points),Time_range),Num_Points,Time_range))
title(main="Second Principal Field")

win.graph()
plot(t(matrix(rep(1:Time_range, Num_Points),Time_range,Num_Points)), Principal_Fields[,1,], col=matrix(rep(rainbow(Num_Points),Time_range),Num_Points,Time_range))
title(main="Third Principal Field")
