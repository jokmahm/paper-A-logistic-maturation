Method2 = function(lp1,lp2,up1,up2,l){
  p2L = lp2-l
  p2H = up2-l
  pikma = c(lp1*p2L, lp1*p2H, up1*p2L, up1*p2H)
  p = min(pikma)
  q = max(pikma)
  ymin = 1/(1+exp(4*q))
  ymax = 1/(1+exp(4*p))
  y = c(ymin, ymax)
  return(y)
}
