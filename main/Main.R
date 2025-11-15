# Periodic year 
library("TMB")
compile("mature.cpp")
dyn.load(dynlib("mature"))

# Data
source("Z:/Rdata/maturationData.R")
# Needed functions
source("Z:/R-program/Accumolated-Multistart/Maturing.R")
source("Z:/R-program/Accumolated-Multistart/Method1.R")
source("Z:/R-program/Accumolated-Multistart/Method2.R")

min_age = 2
max_age = 3

#starting = c(1972, 1989, 1991, 1994, 1998, 2004, 2007)
#ending   = c(1987, 1991, 1994, 1998, 2004, 2007, 2011) 

#starting = 1972:2018
#ending   = 1973:2019

starting = 1989
ending   = 1990

ny <- length(starting)

p1 <- matrix(0,ny,2)
p2 <- matrix(0,ny,2)
p3 <- matrix(0,ny,2)
nu <- matrix(0,ny,2)
y  <- matrix(0,ny,6)

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
  for (r in seq(0.1,1,0.1)){
    for (s in seq(12,16,0.5)){
      for (m in seq(0.01,0.1,0.01)){
        for (n in seq(1,100,1)){
          parameters <- list(lnp1=log(r),
                             lnp2=log(s),
                             lnp3=log(m),
                             lnnu=log(n))
          obj <- MakeADFun(data, parameters, DLL = "mature")
          opt <- nlminb(obj$par, obj$fn, obj$gr
                  #,control = list(eval.max = 1000, iter.max = 500)
                  #,lower = c(log(0.01), log(12), log(10^(-15)), -Inf) 
                  #,upper = c(log(5), log(16), log(0.2), Inf)
          )

          rep <- sdreport(obj)
          Epar <- summary(rep)

          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.06 & Epar[6,1]>12 & Epar[6,1]<16 & num==1){
            y[i,1] = start_year
            y[i,2] = r
            y[i,3] = s
            y[i,4] = m
            y[i,5] = n
            y[i,6] = opt$convergence
            p1[i,] = Epar[5,]
            p2[i,] = Epar[6,]
            p3[i,] = Epar[7,]
            nu[i,] = Epar[8,]
            num = 2
            }
          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.06 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){
            y[i,1] = start_year
            y[i,2] = r
            y[i,3] = s
            y[i,4] = m
            y[i,5] = n
            y[i,6] = opt$convergence
            p1[i,] = Epar[5,]
            p2[i,] = Epar[6,]
            p3[i,] = Epar[7,]
            nu[i,] = Epar[8,]
            break}
        }
        if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.06  & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
      }
      if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.06  & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
    }
    if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE  & Epar[5,1]>0.06  & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
  }
}
Epar
# Correct years that optimization didn't converge
if (min_age == 2){
  p1[45,1] = (p1[44,1]+p1[46,1])/2
  p2[45,1] = (p2[44,1]+p2[46,1])/2
  p3[45,1] = (p3[44,1]+p3[46,1])/2} # Age 2-3
if (min_age == 3){
  p1[16,1] = (p1[15,1]+p1[17,1])/2
  p2[16,1] = (p2[15,1]+p2[17,1])/2
  p3[16,1] = (p3[15,1]+p3[17,1])/2
  #p1[26,] = c(0.02010468, 0.0004329296)
  #p2[26,] = c(14.19594, 0.02154666)
  #p3[26,] = c(0.02482049, 0.0001021981)
  #nu[26,] = c(249903280631, 5658370642)
  #y[26,] = c(1997, 0.1, 12, 0.01, 5, 1)
  #p1[36,] = c(0.02175657, 4.501866e-03)
  #p2[36,] = c(15.72158, 1.361757e-01)
  #p3[36,] = c(0.002097124, 1.729941e-03)
  #nu[36,] = c(1.718303e+09, 1866752304)
  #y[36,] = c(2007, 0.1, 15.0, 0.03,  6, 1)
  }  # Age 3-4

# plotting the results
plot(starting,p1[,1],type='b', pch = 19, col = "red", xlab = "Years", ylab = "P1",main = 'variations')
plot(starting,p2[,1],type='b', pch = 19, col = "red", xlab = "Years", ylab = "P2",main = 'variations')
plot(starting,p3[,1],type='b', pch = 19, col = "red", xlab = "Years", ylab = "P3",main = 'variations')
#____________________________________________________________________________
#____________________________________________________________________________
# maturation curves

length_l = D$length_l
meanlength = D$meanlength

r = Maturing(meanlength,p1[1,1],p2[1,1])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1))
if (ny>1){
  for(i in 2:ny){
    r = Maturing(meanlength,p1[i,1],p2[i,1])
    lines(meanlength,r,type = 'l',col = i)
  }
  legend("topleft", legend = paste(starting,'-',ending), col=1:ny, pch='__', cex=0.7)
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
TB = TB/max(TB)

plotyy(years,TB,years,p1[,1], lty = 1, xlab = "Year", ylab = "Total Biomass", main = "P1 & Total Biomass")
mtext("P1",line = 1, side=4)
plotyy(years,TB,years,p2[,1], lty = 1, xlab = "Year", ylab = "Total Biomass", main = "P2 & Total Biomass")
mtext("P2",line = 1, side=4)
plotyy(years,TB,years,p3[,1], lty = 1, xlab = "Year", ylab = "Total Biomass", main = "P3 & Total Biomass")
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
k = 5   # frame size for moving average
ma1 = MA(p1[,1],k)
ma2 = MA(p2[,1],k)
ma3 = MA(p3[,1],k)
plotyy(years,TB,years,ma1, lty = 5, xlab = "Year", ylab = "Total Biomass", main = "moving average P1 & Total Biomass")
mtext("moving average P1",line = 1, side=4)
plotyy(years,TB,years,ma2, lty = 5, xlab = "Year", ylab = "Total Biomass", main = "moving average P2 & Total Biomass")
mtext("moving average P2",line = 1, side=4)
plotyy(years,TB,years,ma3, lty = 5, xlab = "Year", ylab = "Total Biomass", main = "moving average P3 & Total Biomass")
mtext("moving average P3",line = 1, side=4)

#____________________________________________________________________________
# Classification of the maturation functions 

# group 1 'horizontal'
if (min_age == 2){c = c(15,18,39,40)}    # Age 2-3
if (min_age == 3){c = c(14,15,21,22,24,33,44)}  # Age 3-4
c1 = p1[c,1]
c2 = p2[c,1]
r = Maturing(meanlength,c1[1],c2[1])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 1')
for(i in 2:length(c)){
  r = Maturing(meanlength,c1[i],c2[i])
  lines(meanlength,r,type = 'l',col = i)
}
legend("topleft", legend = paste(starting[c],'-',ending[c]), col=1:length(c), pch='__', cex=0.7)

# group 2 'vertical'
if (min_age == 2){d = c(7,13,14,19,22,23,31,35)}    # Age 2-3
if (min_age == 3){d = c(4,16,19,26,36,39)}  # Age 3-4
d1 = p1[d,1]
d2 = p2[d,1]
r = Maturing(meanlength,d1[1],d2[1])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 2')
for(i in 2:length(d)){
  r = Maturing(meanlength,d1[i],d2[i])
  lines(meanlength,r,type = 'l',col = i)
}
legend("topleft", legend = paste(starting[d],'-',ending[d]), col=1:length(d), pch='__', cex=0.7)

# group 3 'normal'
e = 1:length(p1[,1])
e = setdiff(e,c)
e = setdiff(e,d)
e1 = p1[e,1]
e2 = p2[e,1]
r = Maturing(meanlength,e1[1],e2[1])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 3')
for(i in 2:length(e)){
  r = Maturing(meanlength,e1[i],e2[i])
  lines(meanlength,r,type = 'l',col = i)
}
legend("topleft", legend = paste(starting[e],'-',ending[e]), col=1:length(d), pch='__', cex=0.7)

u = rep(0,47)
u[c] = 1
u[d] = 2
u[e] = 3
dat <- data.frame(x=years, y1 = TB)
plot(y1 ~ x, data=dat, type = 'n', xlab = 'Year', ylab = 'Total Biomass', main = 'Age 2-3', ylim = c(0, 1))     
text(dat$x,dat$y1,label = u, col=4)
lines(years, TB, type = 'l', col = 3)

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
clusplot(df34,kmf34$cluster,main = 'Clustring maturation results', xlab = 'P1', ylab = 'P2')


for(i in 1:length(years)){
  if (result34[i,2]==1){
    cont = i
    break
    }
}
r = Maturing(meanlength,result34[cont,3],result34[cont,4])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 1')
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
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 2')
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
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 3')
for(i in 2:length(years)){
  if (result34[i,2]==3){
    r = Maturing(meanlength,result34[i,3],result34[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}

dat <- data.frame(x=years, y1=TB)
plot(y1 ~ x, data=dat, type='n', xlab = 'Year', ylab = 'Total Biomass', main = 'Age 3-4', ylim=c(0, 1))     
text(dat$x,dat$y1,label=result34[,2],col=1)
lines(years,TB,type = 'l',col = 6)

#_____________________________________________________
# Clustering Age-23

dat23 = cbind(p1[,1],p2[,1])
df23 = scale(dat23)
kmf23 = kmeans(df23,3)
res23 = cbind(kmf23$cluster,dat23)
result23 = cbind(1973:2019,res23)
clusplot(df23,kmf23$cluster,main = 'Clustring maturation results', xlab = 'P1', ylab = 'P2')


for(i in 1:length(years)){
  if (result23[i,2]==1){
    cont = i
    break
  }
}
r = Maturing(meanlength,result23[cont,3],result23[cont,4])
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 1')
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
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 2')
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
plot(meanlength,r,type = 'l',col = 1,xlab='mean length',ylab='maturation',ylim=c(0,1), main = 'group 3')
for(i in 2:length(years)){
  if (result23[i,2]==3){
    r = Maturing(meanlength,result23[i,3],result23[i,4])
    lines(meanlength,r,type = 'l',col = i)
  }  
}


dat <- data.frame(x=years, y1=TB)
plot(y1 ~ x, data=dat, type='n', xlab = 'Year', ylab = 'Total Biomass', main = 'Age 2-3', ylim=c(0, 1))     
text(dat$x,dat$y1,label=result23[,2],col=1)
lines(years,TB,type = 'l',col = 3)

final=cbind(result34,res23)
#______________________________________________________