update_sigma2 <- function(lambda, df, Ak_draws, Mu_cond) {

  num_total <- 0
  den_count <- 0

  units <- sort(unique(df$unit))

  for (k in seq_along(units)) {

    dfk <- df[df$unit == units[k], ]
    t <- dfk$t
    Y <- dfk$Y
    n <- length(Y)
    dlt <- diff(t)

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

    num <- (Y[-1] - mu_E)^2
    den <- 1 - exp(-2 * lambda * dlt)

    num_total <- num_total + sum(num / den)
    den_count <- den_count + length(num)
  }

  sigma2_new <- (2 * lambda / den_count) * num_total
  sigma2_new
}
