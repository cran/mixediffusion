posterior_ou <- function(ak, dfk, theta, Mu_cond){

  t     <- dfk$t
  #dlt_t <- diff(t)
  #dlt_Y <- diff(dfk$Y)

  mu_a     <- theta["mu_a"]
  sigma2_a <- theta["sigma2_a"]
  sigma2   <- theta["sigma2"]
  lambda <- theta["lambda"]
  # --- log prior (hasta constante)
  log_prior <- -0.5 * (
    (ak - mu_a)^2 / sigma2_a +
      log(sigma2_a)
  )

  # --- log likelihood (incrementos)
  mu <- Mu_cond(
    ak   = ak,
    Yi_1 = dfk$Y[-length(t)],
    ti   = t[-1],
    ti_1 = t[-length(t)],
    lambda = lambda
  )
  sigma_cond <- S2cond(S2 = sigma2, lambda = lambda,
                       ti   = t[-1],
                       ti_1 = t[-length(t)])
  log_like <- -0.5 * sum(
    log(sigma_cond) +
      (dfk$Y[-1] - mu)^2 / (sigma_cond)
  )

  return(log_prior + log_like)
}
