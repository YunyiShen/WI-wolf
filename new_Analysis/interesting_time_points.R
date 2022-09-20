# root of 3rd order derivatives
root_3rd_der <- function(r,p0, invK, t0 = 1980){
  ap  <- invK * p0
  inside <- sqrt(3) * sqrt(ap^6 - 4 * ap^5 + 6*ap^4 - 4 * ap^3 + ap^2 ) 
  outside <- -2 * ap^3 +4 * ap^2 - 2 * ap
  
  x1 <- log((inside+outside)/(ap^3-ap^2))/r
  x2 <- log((outside-inside)/(ap^3-ap^2))/r
  return(c(x1,x2)+t0)
  
}

# root of 2nd order derivatives 
root_2nd_der <- function(r, p0, invK, t0 = 1980){
  t0+(-sum(log(c(p0,invK))))/r
}

