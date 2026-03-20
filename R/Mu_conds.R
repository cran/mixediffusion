Mu_cond_sin <- function(ak, Yi_1, ti, ti_1, lambda){
  term_drift <- Yi_1 * exp(-lambda * (ti - ti_1))
  term1 <- lambda / (lambda^2 + ak^2)
  term2 <- lambda * sin(ak * ti) - ak * cos(ak * ti)
  term3 <- exp(-lambda * (ti - ti_1)) *
    (lambda * sin(ak * ti_1) - ak * cos(ak * ti_1))
  value <- term_drift + term1 * (term2 - term3)
  return(value)
}

Mu_cond_cos <- function(ak,Yi_1,ti,ti_1,lambda){
  term_drift <- Yi_1 * exp(-lambda * (ti-ti_1))
  term1 <- lambda/(lambda^2 + ak^2)
  term2 <- lambda*cos(ak*ti)+ak*sin(ak*ti)
  term3 <- exp(-lambda*(ti-ti_1))*(lambda*cos(ak*ti_1)+ak*sin(ak*ti_1))
  value <- term_drift + term1*(term2 - term3)
  return(value)
}


Mu_cond_at1 <- function(ak, Yi_1, ti, ti_1, lambda){
  term_drift <- Yi_1 * exp(-lambda * (ti - ti_1))
  term2 <- ti - 1 / lambda
  term3 <- exp(-lambda * (ti - ti_1)) * (ti_1 - 1 / lambda)
  value <- term_drift + ak * (term2 - term3)
  return(value)
}

Mu_cond_at2 <- function(ak, Yi_1, ti, ti_1, lambda){
  term_drift <- Yi_1 * exp(-lambda * (ti - ti_1))
  term2 <- ti^2 - 2 * ti / lambda + 2 / lambda^2
  term3 <- exp(-lambda * (ti - ti_1)) *
    (ti_1^2 - 2 * ti_1 / lambda + 2 / lambda^2)
  value <- term_drift + ak * (term2 - term3)
  return(value)
}


Mu_cond_exp <- function(ak, Yi_1, ti, ti_1, lambda){
  term_drift <- Yi_1 * exp(-lambda * (ti - ti_1))
  term1 <- lambda / (lambda + ak)
  term2 <- exp(ak * ti)
  term3 <- exp(-lambda * (ti - ti_1)) * exp(ak * ti_1)
  value <- term_drift + term1 * (term2 - term3)
  return(value)
}

Mu_cond_const <- function(ak, Yi_1, ti, ti_1, lambda){
  delta_t <- ti - ti_1
  term_drift <- Yi_1 * exp(-lambda * delta_t)
  value <- term_drift + ak * (1 - exp(-lambda * delta_t))
  return(value)
}



