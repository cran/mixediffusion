posterior_winner <- function(ak, dfk, theta, Mu_dlt){

  t     <- dfk$t
  dlt_t <- diff(t)
  dlt_Y <- diff(dfk$Y)

  mu_a     <- theta["mu_a"]
  sigma2_a <- theta["sigma2_a"]
  sigma2   <- theta["sigma2"]

  # --- log prior (hasta constante)
  log_prior <- -0.5 * (
    (ak - mu_a)^2 / sigma2_a +
      log(sigma2_a)
  )

  # --- log likelihood (incrementos)
  mu_inc <- Mu_dlt(
    ak   = ak,
    ti   = t[-1],
    ti_1 = t[-length(t)]
  )

  log_like <- -0.5 * sum(
    log(sigma2 * dlt_t) +
      (dlt_Y - mu_inc)^2 / (sigma2 * dlt_t)
  )

  return(log_prior + log_like)
}
