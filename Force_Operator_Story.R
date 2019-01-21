source("C:/Users/tpl97/Documents/R/Force_Operator.R")

Num_Points <- 6;
Time_range <- 100;
delta_Time <- .1

Position <- array(NA, c(Num_Points,Time_range))
Velocity <- array(NA, c(Num_Points,Time_range))
Force_vector <- array(NA, c(Num_Points,Time_range))

Position_Spectrum <- array(NA, c(Num_Points,Time_range))
Velocity_Spectrum <- array(NA, c(Num_Points,Time_range))
Force_Spectrum <- array(NA, c(Num_Points,Time_range))

Principal_Fields <- array(NA, c(Num_Points,Num_Points,Time_range))

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
f = function(r) -r^(-2)*1000 #Gravity
# f = function(r) -1/r*100 #1/r

while(t < Time_range) {
  
  t = t + 1
  
  F <- Force_Operator(z,f)
  Vac <- matrix(rep(1,Num_Points),Num_Points,1)
  Force_vector[,t] <- -F %*% Vac
  
  #Building "Force Hamiltonian" Operator
  Force_Hamiltonian <- complex(real = 0,imaginary = 1)*F

  #Finding "Force Hamiltonian" Spectrum
  Force_eigenstructure <- eigen(Force_Hamiltonian, symmetric=TRUE, only.values = FALSE)
  
  Position[,t] <- z
  Velocity[,t] <- Vz
  
  Principal_Fields[,,t] <- Force_eigenstructure$vectors
  Position_Spectrum[,t] <- diag(Conj(t(Force_eigenstructure$vectors)) %*% diag(Position[,t]) %*% Force_eigenstructure$vectors)
  Velocity_Spectrum[,t] <- diag(Conj(t(Force_eigenstructure$vectors)) %*% diag(Velocity[,t]) %*% Force_eigenstructure$vectors)
  Force_Spectrum[,t] <- Force_eigenstructure$values
  
  #Time Evoloution
  z = z + Vz * delta_Time
  Vz = Vz + Force_vector[,t] * delta_Time
  
}

win.graph()
plot(Position, col=rainbow(Num_Points))
title(main="Position")
win.graph()
plot(Velocity, col=rainbow(Num_Points))
title(main="Velocity")
win.graph()
plot(Force_vector, col=rainbow(Num_Points))
title(main="Force")

win.graph()
plot(Position_Spectrum)
title(main="Position Spectrum")
win.graph()
plot(Velocity_Spectrum)
title(main="Velocity Spectrum")
win.graph()
# For SVD:
plot(t(matrix(rep(1:Time_range, Num_Points),Time_range,Num_Points)), Force_Spectrum, col=matrix(rep(rainbow(Num_Points),Time_range),Num_Points,Time_range))
title(main="Force Spectrum")

if (FALSE) {
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

