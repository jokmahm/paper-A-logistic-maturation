Method1 = function(r,p1,p2,sdp1,sdp2,l){
  vp1 = sdp1^2;
  vp2 = sdp2^2;

  v = ((4*exp(4*p1*(p2-l)))*r^2)^2*((p2-l)^2*vp1+p1^2*vp2);
  return(v)
}
