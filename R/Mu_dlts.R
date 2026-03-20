Mu_dlt_cos <- function(ak,ti,ti_1){
  value <- (1/ak)*(sin(ak*ti)-sin(ak*ti_1))
  return(value)
}

Mu_dlt_exp <- function(ak,ti,ti_1){
  value <- (1/ak)*(exp(ak*ti)-exp(ak*ti_1))
  return(value)
}

Mu_dlt_power <- function(ak, ti, ti_1, p){
  if (abs(p + 1) < 1e-8) {
    # p = -1
    ak * (log(ti) - log(ti_1))
  } else {
    ak / (p + 1) * (ti^(p + 1) - ti_1^(p + 1))
  }
}

Mu_dlt_sin <- function(ak, ti, ti_1){
  value <- (1 / ak) * (cos(ak * ti_1) - cos(ak * ti))
  return(value)
}


