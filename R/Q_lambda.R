Q_lambda <- function(lambda, df, Ak_draws, sigma2, Mu_cond) {

  Q_total <- 0
  units <- sort(unique(df$unit))

  for (k in seq_along(units)) {

    dfk <- df[df$unit == units[k], ]
    t <- dfk$t
    Y <- dfk$Y
    n <- length(Y)
    dlt <- diff(t)

    # matriz (n-1) x M de medias condicionales para cada draw de ak
    mu_mat <- sapply(
      Ak_draws[k, ],
      function(ak_m)
        Mu_cond(
          ak     = ak_m,
          Yi_1  = Y[-n],
          ti    = t[-1],
          ti_1  = t[-n],
          lambda = lambda
        )
    )

    mu_E <- rowMeans(mu_mat)

    sigma_cond <- (sigma2 / (2 * lambda)) * (1 - exp(-2 * lambda * dlt))

    Q_total <- Q_total + (-0.5) * sum(
      log(sigma_cond) + (Y[-1] - mu_E)^2 / sigma_cond
    )
  }

  Q_total
}
