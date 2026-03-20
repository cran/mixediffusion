S2cond <- function(S2,lambda,ti,ti_1){
  dlt <- ti-ti_1
  value <- (S2/(2*lambda))*(1-exp(-2*lambda*dlt))
  return(value)
}
