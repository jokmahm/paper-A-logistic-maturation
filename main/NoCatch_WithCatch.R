# Periodic year 
library("TMB")
compile("mature.cpp")
dyn.load(dynlib("mature"))

# Data
source("Z:/Rdata/maturationData.R")

min_age = 3
max_age = 4

# Years with Catch
#____________________________________________________________________________

starting = 1972
ending   = 1980

ny <- 1

p1 <- matrix(0,ny,2)
p2 <- matrix(0,ny,2)
p3 <- matrix(0,ny,2)
nu <- matrix(0,ny,2)
y  <- matrix(0,ny,6)

for (i in 1:ny){
  num = 1
  start_year = starting[i]
  end_year = ending[i]
  
  D = maturationData(1988,1990,min_age,max_age)
  R = D$Nl
  K = D$N
  E = D$C
  D = maturationData(1994,1997,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  D = maturationData(2004,2006,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  D = maturationData(2016,2017,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  data <- list(min_age = D$min_age,
               max_age = D$max_age,
               start_year = start_year,
               end_year = end_year,
               length_l = D$length_l,
               meanlength = D$meanlength,
               Nl = R,
               N = K,
               C = E)
  for (r in seq(0.1,1,0.1)){
    for (s in seq(12,16,0.5)){
      for (m in seq(0.01,0.1,0.01)){
        for (n in seq(1,10,1)){
          parameters <- list(lnp1=log(r),
                             lnp2=log(s),
                             lnp3=log(m),
                             lnnu=log(n))
          obj <- MakeADFun(data, parameters, DLL = "mature")
          opt <- nlminb(obj$par, obj$fn, obj$gr)

          rep <- sdreport(obj)
          Epar <- summary(rep)
          
          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & num==1){
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
          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){
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
        if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break} #& p1[i,]>0.01 & p1[i,]<1 & p2[i,]>13 & p2[i,]<15
      }
      if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
    }
    if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE  &  Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
  }
}
Epar
E1 = Epar[5:8,]
# plotting the results
#____________________________________________________________________________
lp1 = p1[,1]-p1[,2]
up1 = p1[,1]+p1[,2]
lp2 = p2[,1]-p2[,2]
up2 = p2[,1]+p2[,2]

source("Z:/R-program/Accumolated-Multistart/Maturing.R")
source("Z:/R-program/Accumolated-Multistart/Method1.R")
source("Z:/R-program/Accumolated-Multistart/Method2.R")

meanlength = D$meanlength
r = Maturing(meanlength,p1[1,1],p2[1,1])
plot(meanlength,r,type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
legend("topleft", legend = paste(starting,'-',ending), col=1, pch='__', cex=0.7)


# r1 = Maturing(meanlength,0.2960653648 ,12.0575447019 )

UncerMeth = 2
length_l = D$length_l

if (UncerMeth==1){
  for(i in 1:ny){
    r = Maturing(meanlength,p1[i,1],p2[i,1])
    plot(meanlength,r,type = 'l',xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),main='Years with no catch',ylim=c(0,1))
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
    plot(meanlength,r,type = 'l',xlab='mean length',ylab='maturation',main='Years with no catch',ylim=c(0,1))
    for(j in 1:length_l){
      w=Method2 (lp1[i],lp2[i],up1[i],up2[i],meanlength[j])
      ymin[j] = w[1]
      ymax[j] = w[2]
    }
    lines(meanlength,ymin,type = 'l',lty = 2)
    lines(meanlength,ymax,type = 'l',lty = 2)
  } 
}


# Years with Catch
#____________________________________________________________________________

starting = 1972
ending   = 2006 #2008

ny <- 1

p1 <- matrix(0,ny,2)
p2 <- matrix(0,ny,2)
p3 <- matrix(0,ny,2)
nu <- matrix(0,ny,2)
y  <- matrix(0,ny,6)

for (i in 1:ny){
  num = 1
  start_year = starting[i]
  end_year = ending[i]
  
  D = maturationData(1972,1987,min_age,max_age)
  R = D$Nl
  K = D$N
  E = D$C
  D = maturationData(1990,1994,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  D = maturationData(1997,2004,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  D = maturationData(2006,2014,min_age,max_age)
  R = rbind(R,D$Nl)
  K = rbind(K,D$N)
  E = rbind(E,D$C)
  data <- list(min_age = D$min_age,
               max_age = D$max_age,
               start_year = start_year,
               end_year = end_year,
               length_l = D$length_l,
               meanlength = D$meanlength,
               Nl = R,
               N = K,
               C = E)
  for (r in seq(0.1,1,0.1)){
    for (s in seq(12,16,0.5)){
      for (m in seq(0.01,0.1,0.01)){
        for (n in seq(1,10,1)){
          parameters <- list(lnp1=log(r),
                             lnp2=log(s),
                             lnp3=log(m),
                             lnnu=log(n))
          obj <- MakeADFun(data, parameters, DLL = "mature")
          opt <- nlminb(obj$par, obj$fn, obj$gr)
          rep <- sdreport(obj)
          Epar <- summary(rep)
          
          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & num==1){
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
          if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){
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
        if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break} #& p1[i,]>0.01 & p1[i,]<1 & p2[i,]>13 & p2[i,]<15
      }
      if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
    }
    if (is.nan(Epar[5,2])==FALSE & is.nan(Epar[6,2])==FALSE & is.nan(Epar[7,2])==FALSE  & Epar[5,1]>0.1 & Epar[6,1]>12 & Epar[6,1]<16 & opt$convergence==0){break}
  }
}
Epar
E2 = Epar[5:8,]
# plotting the results
#____________________________________________________________________________
lp1 = p1[,1]-p1[,2]
up1 = p1[,1]+p1[,2]
lp2 = p2[,1]-p2[,2]
up2 = p2[,1]+p2[,2]

source("Z:/R-program/Accumolated-Multistart/Maturing.R")
source("Z:/R-program/Accumolated-Multistart/Method1.R")
source("Z:/R-program/Accumolated-Multistart/Method2.R")

meanlength = D$meanlength
rnc = Maturing(meanlength,E1[1,1],E1[2,1])
plot(meanlength,rnc, bty='l' ,lwd=3, type = 'l',col = 1,xlab=expression(mu[l]),ylab=expression('r(l,p'[1]*',p'[2]*')'),ylim=c(0,1))
rwc = Maturing(meanlength,E2[1,1],E2[2,1])
lines(meanlength,rwc,lwd=3, type = 'l',col = 2)
legend("topleft", legend = c(expression('r'[a]),expression('r'[b])), lwd=3,col=1:2, pch='__', cex=1)


UncerMeth = 2
length_l = D$length_l

if (UncerMeth==1){
  for(i in 1:ny){
    r = Maturing(meanlength,E2[1,1],E2[2,1])
    plot(meanlength,r, bty='l' ,type = 'l',xlab='mean length',ylab='maturation',main='Years with Fishing',ylim=c(0,1))
    v = Method1(r,E2[1,1],E2[2,1],E2[1,2],E2[2,2],meanlength)
    lines(meanlength,r+sqrt(v),type = 'l',lty = 2)
    lines(meanlength,r-sqrt(v),type = 'l',lty = 2)
  }
}

if (UncerMeth==2){
  for(i in 1:ny){
    ymin = rep(0,length_l);
    ymax = rep(0,length_l);
    r = Maturing(meanlength,E2[1,1],E2[2,1])
    plot(meanlength,r, bty='l' ,type = 'l',xlab='mean length',ylab='maturation',main='Years with Fishing',ylim=c(0,1))
    for(j in 1:length_l){
      w=Method2 (lp1[i],lp2[i],up1[i],up2[i],meanlength[j])
      ymin[j] = w[1]
      ymax[j] = w[2]
    }
    lines(meanlength,ymin,type = 'l',lty = 2)
    lines(meanlength,ymax,type = 'l',lty = 2)
  } 
}

# Age 3-4 estimated alpha and beta

loges = function(a,b,ml){1/(1+exp(a*(ml-b)))}
ap = loges(-1.531304, 13.47015, meanlength)

plot(meanlength,rnc/rwc, bty='l', type = 'l', lwd=3,col = 4,xlab=expression(mu[l]),ylab='Ratio',ylim=c(0,1))
lines(meanlength,ap,lwd=3, lty = 2,col = 2)
legend("topleft", legend = c(expression('r'[a]*'/r'[b]),expression('r'[ab])), lwd=3, col=c(4,2), lty = c(1,2), cex=1)


#plot(meanlength,(rwc-rnc)/meanlength,type = 'l')

#___________________________________________________
library(xtable)
tab1 <- xtable(format(E1, scientific = TRUE, digits=3), caption= "summary statistics of air pollution data")
tab2 <- xtable(format(E2, scientific = TRUE, digits=3), caption= "summary statistics of air pollution data")
