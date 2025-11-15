# Parameter estimator for mature Capelin
library("TMB")
compile("mature.cpp")
dyn.load(dynlib("mature"))

# The estimation parameter will be from 1972 to ''end_simulat_year'
end_simulat_year = 1986

# Data
min_age=3
max_age=4
start_year = 1972

ny <- end_simulat_year-1979
p1 <- matrix(0,ny,2)
p2 <- matrix(0,ny,2)
p3 <- matrix(0,ny,2)
for (i in 1:ny){
  
  end_year = 1979+i
  D <- maturationData(start_year,end_year,min_age,max_age)
  meanlength <- D$meanlength
  length_l <- D$length_l
  Nl <- D$Nl
  N <- D$N
  C <- D$C

  data <- list(min_age = min_age,
                   max_age = max_age,
                   start_year = start_year,
                   end_year = end_year,
                   length_l = length_l,
                   meanlength = meanlength,
                   Nl = Nl,
                   N = N,
                   C = C)
  parameters <- list(lnp1=log(.01),#1.2,
                   lnp2=log(13),#2.6,
                   lnp3=log(.09),#-3.5,
                   lnnu=log(1.3))#2)

  obj <- MakeADFun(data, parameters, DLL = "mature")
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  # Get standard deviation of parameters and latent process
  rep <- sdreport(obj)

  # Get Report
  Epar <- summary(rep)
  p1[i,] = Epar[5,]
  p2[i,] = Epar[6,]
  p3[i,] = Epar[7,]
}

# plotting the results
library(ggplot2)
#____________________________________________________________________________
data.estimates = data.frame(
  var   = 1:ny,
  par = p1[,1],
  se = p1[,2])
data.estimates$idr <- data.estimates$par
data.estimates$upper <- (data.estimates$par + data.estimates$se)
data.estimates$lower <- (data.estimates$par - data.estimates$se)

f1 <- ggplot(data.estimates, aes(var,idr)) 
f1 + geom_point() + 
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years starting 1980") + ylab("P1")  + coord_flip()

#____________________________________________________________________________
data.estimates = data.frame(
  var   = 1:ny,
  par = p2[,1],
  se = p2[,2])
data.estimates$idr <- data.estimates$par
data.estimates$upper <- (data.estimates$par + data.estimates$se)
data.estimates$lower <- (data.estimates$par - data.estimates$se)

f2 <- ggplot(data.estimates, aes(var,idr)) 
f2 + geom_point() + 
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years starting 1980") + ylab("P2") + coord_flip()

#____________________________________________________________________________
data.estimates = data.frame(
  var   = 1:ny,
  par = p3[,1],
  se = p3[,2])
data.estimates$idr <- data.estimates$par
data.estimates$upper <- (data.estimates$par + data.estimates$se)
data.estimates$lower <- (data.estimates$par - data.estimates$se)

f3 <- ggplot(data.estimates, aes(var,idr)) 
f3 + geom_point() +
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years starting 1980") + ylab("P3")  + coord_flip()
#____________________________________________________________________________