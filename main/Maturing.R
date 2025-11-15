Maturing=function(meanlength,p1,p2){
  
  r = 1.0/(1.0+exp(4*p1*(p2-meanlength)))
  
  return(r)
}
