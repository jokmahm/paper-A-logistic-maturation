# Periodic year 
library("TMB")
compile("mature.cpp")
dyn.load(dynlib("mature"))

set.seed(2020)
# set.seed(11)

# Data
codbiomass <- scan("CodBiomass.txt")
codbiomass <- codbiomass/10^6
source("Z:/Rdata/maturationData.R")
# Needed functions
source("Z:/R-program/Accumolated-Multistart/Maturing.R")
source("Z:/R-program/Accumolated-Multistart/Method1.R")
source("Z:/R-program/Accumolated-Multistart/Method2.R")

min_age = 3
max_age = 4

starting = 1972:2018
ending   = 1973:2019

if (min_age == 2){
  y <- read.table("InitialParameterValue23.txt", header = FALSE)
}
if (min_age == 3){
  y <- read.table("InitialParameterValue34.txt", header = FALSE)
}

ny <- length(starting)

p1 <- matrix(0,ny,2)
p2 <- matrix(0,ny,2)
p3 <- matrix(0,ny,2)
nu <- matrix(0,ny,2)

for (i in 1:ny){
  num = 1
  start_year = starting[i]
  end_year = ending[i]
  
  D = maturationData(start_year,end_year,min_age,max_age)
  
  data <- list(min_age = D$min_age,
               max_age = D$max_age,
               start_year = D$start_year,
               end_year = D$end_year,
               length_l = D$length_l,
               meanlength = D$meanlength,
               Nl = D$Nl,
               N = D$N,
               C = D$C)

  parameters <- list(lnp1=log(y[i,2]),
                    lnp2=log(y[i,3]),
                    lnp3=log(y[i,4]),
                    lnnu=log(y[i,5]))
          obj <- MakeADFun(data, parameters, DLL = "mature")
          opt <- nlminb(obj$par, obj$fn, obj$gr)
          
  rep <- sdreport(obj)
  Epar <- summary(rep)
          
  p1[i,] = Epar[5,]
  p2[i,] = Epar[6,]
  p3[i,] = Epar[7,]
  nu[i,] = Epar[8,]
}
          
# Correct years that optimization didn't converge
if (min_age == 2){
  p1[45,] = (p1[44,]+p1[46,])/2
  p2[45,] = (p2[44,]+p2[46,])/2
  p3[45,] = (p3[44,]+p3[46,])/2} # Age 2-3
if (min_age == 3){
  p1[16,] = (p1[15,]+p1[17,])/2
  p2[16,] = (p2[15,]+p2[17,])/2
  p3[16,] = (p3[15,]+p3[17,])/2
  
  p1[36,] = (p1[35,]+p1[37,])/2
  p2[36,] = (p2[35,]+p2[37,])/2
  p3[36,] = (p3[35,]+p3[37,])/2}  # Age 3-4

# plotting the estimated parameters

plot(starting,p1[,1],type='b', pch = 19, bty="l", col = "red", xlab = "Year", ylab = expression(paste("Maturation intensity (p"[1],")")))
plot(starting,p2[,1],type='b', pch = 19, bty="l", col = "red", xlab = "Year", ylab = expression(paste("Median length at maturity (p"[2],")")))
plot(starting,p3[,1],type='b', pch = 19, bty="l", col = "red", xlab = "Year", ylab = expression(paste("Mortality (p"[3],")")))

# 2 axis plot (p1-p2 in one graph)
st <- seq(0,46,1)
y <- list(p1[,1], p2[,1])
par(oma = c(0, 2, 2, 3))
plot(st, y[[1]], yaxt = "n", xaxt = "n", lwd = 2, type='l', col = "red", bty="l", xlab = "Year", ylab = "")
axis(1, seq(0,46,1),las=2)
lines(st, y[[1]], lwd = 2, col ="red")
axis(at = pretty(y[[1]]), col.axis = 'red', cex.axis = 1.1, side = 2, col ="red")
sides <- list(4) 
par(new = TRUE)
plot(st, y[[2]], axes = FALSE, type='l', bty="l", col = "blue", xlab = "", ylab = "")
axis(at = pretty(y[[2]]), side = 4, col.axis = 'blue', cex.axis = 1.1, line = NA, col = "blue")
lines(st, y[[2]], lwd = 2, type='l', col = "blue")
mtext(expression("p"[1]),line = 2, side=2, cex = 1.1, col = 'red')
mtext(expression("p"[2]),line = 2, side=4, cex = 1.1, col = 'blue')

plot(st,p3[,1],type='l', pch = 19, lwd = 2,  cex.lab = 1.1, xaxt = "n", bty="l", col = "steelblue", xlab = "Year", ylab = expression("p"[3]))
axis(1, seq(0,46,1),las=2)

#____________________________________________________________________________
#____________________________________________________________________________
# maturation curves

length_l = D$length_l
meanlength = D$meanlength

r = Maturing(meanlength,p1[1,1],p2[1,1])
plot(meanlength,r,type = 'l',col = 1,xlab='Length',ylab='Maturation',ylim=c(0,1))
if (ny>1){
  for(i in 2:ny){
    r = Maturing(meanlength,p1[i,1],p2[i,1])
    lines(meanlength,r,type = 'l',col = i)
  }
  #legend("topleft", legend = paste(starting,'-',ending), col=1:ny, pch='__', cex=0.7)
}

#____________________________________________________________________________
# Uncertainty plot

UncerMeth = 2           # select a method 1:=variance sensitivity  & 2:=Interval calculation

lp1 = p1[,1]-p1[,2]
up1 = p1[,1]+p1[,2]
lp2 = p2[,1]-p2[,2]
up2 = p2[,1]+p2[,2]

if (UncerMeth==1){
  for(i in 1:ny){
    r = Maturing(meanlength,p1[i,1],p2[i,1])
    plot(meanlength,r,type = 'l',xlab='mean length',ylab='maturation',main=paste(starting[i],'-',ending[i]),ylim=c(0,1))
    v = Method1(r,p1[i,1],p2[i,1],p1[i,2],p2[i,2],meanlength)
    lines(meanlength,r+sqrt(v),type = 'l',lty = 2)
    lines(meanlength,r-sqrt(v),type = 'l',lty = 2)
  }
}

if (UncerMeth==2){
  for(i in 1:ny){
    ymin = rep(0,length_l);
    ymax = rep(0,length_l);
    r = Maturing(meanlength,p1[i,1],p2[i,1])
    plot(meanlength,r,type = 'l',xlab='mean length',ylab='maturation',main=paste(starting[i],'-',ending[i]),ylim=c(0,1))
    for(j in 1:length_l){
      w=Method2 (lp1[i],lp2[i],up1[i],up2[i],meanlength[j])
      ymin[j] = w[1]
      ymax[j] = w[2]
    }
    lines(meanlength,ymin,type = 'l',lty = 2)
    lines(meanlength,ymax,type = 'l',lty = 2)
  } 
}

#____________________________________________________________________________
#____________________________________________________________________________
# ggplot parameters and estimated Standard Deviations

library(ggplot2)
# p1

data.estimates = data.frame(
  var   = starting,
  par = p1[,1],
  se = p1[,2])
data.estimates$idr <- data.estimates$par
data.estimates$upper <- (data.estimates$par + data.estimates$se)
data.estimates$lower <- (data.estimates$par - data.estimates$se)

f1 <- ggplot(data.estimates, aes(var,idr)) 
f1 + geom_point() + 
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years") + ylab("P1")  #+ coord_flip()

#____________________________________________________________________________
# p2

data.estimates = data.frame(
  var   = starting,
  par = p2[,1],
  se = p2[,2])
data.estimates$idr <- data.estimates$par
data.estimates$upper <- (data.estimates$par + data.estimates$se)
data.estimates$lower <- (data.estimates$par - data.estimates$se)

f2 <- ggplot(data.estimates, aes(var,idr)) 
f2 + geom_point() + 
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years") + ylab("P2") #+ coord_flip()

#____________________________________________________________________________
# p3

data.estimates = data.frame(
  var   = starting,
  par = p3[,1],
  se = p3[,2])
data.estimates$idr <- data.estimates$par
ub <- (data.estimates$par + data.estimates$se)
ub[ub>1] <- 1
data.estimates$upper <- ub
lb <- (data.estimates$par - data.estimates$se)
lb[lb<0] <- 0
data.estimates$lower <- lb

f3 <- ggplot(data.estimates, aes(var,idr)) 
f3 + geom_point() +
  geom_errorbar(aes(x = var, ymin = lower, ymax = upper), width = 0.2) + 
  theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  xlab("Years") + ylab("P3")  #+ coord_flip()

#____________________________________________________________________________
#____________________________________________________________________________
# Biomass against parameters (plotyy)

library(plotly)
library(pracma)
source("Z:/Rdata/BiomassData.R")
years = 1973:2019
G = BiomassData(1973,2019)
TB = G$B
TB = TB/1000     # 10^6 t

plotyy(years,TB,years,p1[,1], lty = 4, xlab = "Year", ylab = "Capelin Biomass", main = "P1 & Total Biomass")
mtext("P1",line = 1, side=4)
plotyy(years,TB,years,p2[,1], lty = 4, xlab = "Year", ylab = "Capelin Biomass", main = "P2 & Total Biomass")
mtext("P2",line = 1, side=4)
plotyy(years,TB,years,p3[,1], lty = 4, xlab = "Year", ylab = "Capelin Biomass", main = "P3 & Total Biomass")
mtext("P3",line = 1, side=4)

plotyy(years,codbiomass,years,p3[,1], lty = 2, xlab = "Year", ylab = "Cod Biomass", main = "P3 & Cod Biomass")
mtext("P3",line = 1, side=4)

#_____________________________________________________________________________________
# Moving Average (plotyy)

MA <- function(x, n = 5){
  g = stats::filter(x, rep(1 / n, n), sides = 2)
  m = length(x)
  for (i in 1:(n-1)){g[i] = x[i]}
  for (i in (m-n-1):m){g[i] = x[i]}
  return(g)
}
k = 3   # frame size for moving average
ma1 = MA(p1[,1],k)
ma2 = MA(p2[,1],k)
ma3 = MA(p3[,1],k)
plotyy(years,TB,years,ma1, lwd = 2, lty = 1, gridp = F, xlab = "Year", ylab = expression(paste("Capelin biomass (10"^6,"t)")))
mtext(expression(paste("p"[1]," moving average")),line = 1, side=4)
plotyy(years,TB,years,ma2, lwd = 2, lty = 1, gridp = F, xlab = "Year", ylab = expression(paste("Capelin biomass (10"^6,"t)")))
mtext(expression(paste("p"[2]," moving average")),line = 1, side=4)
plotyy(years,TB,years,ma3, lwd = 2, lty = 1, gridp = F, xlab = "Year", ylab = expression(paste("Capelin biomass (10"^6,"t)")))
mtext(expression(paste("p"[3]," moving average")),line = 1, side=4)
 
plotyy(years,codbiomass,years,ma3, lty = 1, xlab = "Year", ylab = expression(paste("Cod biomass (10"^6,"t)")))
mtext(expression(paste("p"[3]," moving average")),line = 1, side=4)
#_________________________________________________________
# K-mean clustering

library(rattle)
library(cluster)

#_____________________________________________________
# Clustering Age-34
years = 1973:2019

dat34 = cbind(p1[,1],p2[,1])
df34 = scale(dat34)
kmf34 = kmeans(df34,3)
res34 = cbind(kmf34$cluster,dat34)
result34 = cbind(years,res34)
clusplot(df34,kmf34$cluster, bty='l', lwd=1, cex = 1.2, color= TRUE, labels= 4, main = NULL, sub =NULL,  xlab = expression("p"[1]), ylab = expression("p"[2]))


for(i in 1:length(years)){
  if (result34[i,2]==1){
    cont = i
    break
  }
}
r = Maturing(meanlength,result34[cont,3],result34[cont,4])
plot(meanlength,r,bty='l', bty='l', type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result34[i,2]==1){
    r = Maturing(meanlength,result34[i,3],result34[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

for(i in 1:length(years)){
  if (result34[i,2]==2){
    cont = i
    break
  }
}
r = Maturing(meanlength,result34[cont,3],result34[cont,4])
plot(meanlength,r,bty='l', bty='l' , type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result34[i,2]==2){
    r = Maturing(meanlength,result34[i,3],result34[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

for(i in 1:length(years)){
  if (result34[i,2]==3){
    cont = i
    break
  }
}
r = Maturing(meanlength,result34[cont,3],result34[cont,4])
plot(meanlength,r, bty='l' , type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result34[i,2]==3){
    r = Maturing(meanlength,result34[i,3],result34[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

dat <- data.frame(x=years, y1=TB)
plot(y1 ~ x, data=dat, bty='l', type='n', xlab = 'Year', ylab = expression(paste("Capelin biomass (10"^6,"t)")), ylim=c(0, 8))     
text(dat$x,dat$y1,label=result34[,2],col=1)
lines(years,TB, type = 'l',col = 6)


plotyy(years, TB, years,codbiomass, lty = 1, xlab = 'Year', ylab = expression(paste("Capelin Biomass (10"^6,"t)")))
text(dat$x,dat$y1,label = result34[,2], col = 1)
mtext(expression(paste("Cod Biomass (10"^6,"t)")),line = 1, side=4)

dat <- data.frame(x=years, y1=codbiomass)
plot(y1 ~ x, data=dat, type='n', xlab = 'Year', ylab = expression(paste("Cod biomass (10"^6,"t)")))     
text(dat$x,dat$y1,label=result34[,2],col=1)
lines(years,codbiomass,type = 'l',col = 5)

#_____________________________________________________
# Clustering Age-23

dat23 = cbind(p1[,1],p2[,1])
df23 = scale(dat23)
kmf23 = kmeans(df23,3)
res23 = cbind(kmf23$cluster,dat23)
result23 = cbind(1973:2019,res23)
clusplot(df23,kmf23$cluster,bty='l', lwd=1, cex = 1.2, color= TRUE, labels= 4, main = NULL, sub =NULL,  xlab = expression("p"[1]), ylab = expression("p"[2]))


for(i in 1:length(years)){
  if (result23[i,2]==1){
    cont = i
    break
  }
}
r = Maturing(meanlength,result23[cont,3],result23[cont,4])
plot(meanlength,r, bty='l', type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result23[i,2]==1){
    r = Maturing(meanlength,result23[i,3],result23[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

for(i in 1:length(years)){
  if (result23[i,2]==2){
    cont = i
    break
  }
}
r = Maturing(meanlength,result23[cont,3],result23[cont,4])
plot(meanlength,r,bty='l',type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result23[i,2]==2){
    r = Maturing(meanlength,result23[i,3],result23[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

for(i in 1:length(years)){
  if (result23[i,2]==3){
    cont = i
    break
  }
}
r = Maturing(meanlength,result23[cont,3],result23[cont,4])
plot(meanlength,r,bty='l',type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:length(years)){
  if (result23[i,2]==3){
    r = Maturing(meanlength,result23[i,3],result23[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}


dat <- data.frame(x=years, y1=TB)
plot(y1 ~ x, data=dat, bty='l', type='n', xlab = 'Year',ylab = expression(paste("Capelin biomass (10"^6,"t)")), ylim=c(0, 8))     
text(dat$x,dat$y1,label=result23[,2],col=1)
lines(years,TB,type = 'l',col = 3)

final=cbind(result34,res23)
#______________________________________________________
# Markov Chain Analysis

library(markovchain)
library(dplyr)
# age 3-4
verifyMarkovProperty(result34[,2])
# age 2-3
verifyMarkovProperty(result23[,2])

#________________________________________
# Extract Latex table of the Parameters

library(xtable)
pprint <- cbind(p1,cbind(p2,p3))
tab<-xtable(format(pprint, scientific = TRUE, digits=3), caption= "summary statistics of air pollution data")


#________________________________________
# Visualization of maturation intensity

pintensity = c(0.06,0.12,0.24,0.48,0.96)
r = Maturing(meanlength,pintensity[1],14)
plot(meanlength,r,lwd = 3,type = 'l',bty='l', col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:5){
  r = Maturing(meanlength,pintensity[i],14)
  lines(meanlength,r,lwd = 3,type = 'l',col = i)
}
leg <- c( expression('p'[1]*'=0.06'),
          expression('p'[1]*'=0.12'),
          expression('p'[1]*'=0.24'),
          expression('p'[1]*'=0.48'),
          expression('p'[1]*'=0.96')
       )
legend("topleft", legend = leg, lwd = 3, col=1:5, pch='__', cex=1)

# Visualization of median length of maturation

pmedian = c(12,13,14,15,16)
r = Maturing(meanlength,0.24,pmedian[1])
plot(meanlength,r,lwd = 3, type = 'l',bty='l', col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
for(i in 2:5){
  r = Maturing(meanlength,0.24,pmedian[i])
  lines(meanlength,r, lwd = 3,type = 'l',col = i)
}
leg <- c( expression('p'[2]*'=12'),
          expression('p'[2]*'=13'),
          expression('p'[2]*'=14'),
          expression('p'[2]*'=15'),
          expression('p'[2]*'=16')
)
legend("topleft", legend = leg, lwd = 3, col=1:5, pch='__', cex=1)

