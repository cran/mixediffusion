fit_wiener <- function(df,mu = "at^1",tol = 1e-4,
                       max_iter = 100,theta = NULL,
                       M = 100,verbose = TRUE,
                       mu_dlt = NULL,n_mcmc = 1000, burnin = 500){
  if(!is.null(mu_dlt)){
    Mu_dlt <- mu_dlt
  }else{
    Mu_dlt <- make_Mu_dlt(mu)
  }
  log_msg <- function(...) {
    if (isTRUE(verbose)) {
      cat(...)
    }
  }

  if(is.null(theta)){
    theta <- c(
      mu_a     = 1,
      sigma2_a = 1,
      sigma2   = 1
    )
  }
  theta_hist <- list()
  theta_hist[[1]] <- theta
  K <- length(unique(df$unit))
  count_em <- 1
  error <- 1e300

  while (error > tol & count_em < max_iter){
    Ak <- matrix(ncol = M, nrow = K)
    counter <- 1

    num_sigma <- 0
    den_sigma <- 0

    for(k in unique(df$unit)){
      dfk <- df %>%
        filter(unit == k)

      log_post_ak <- function(a){
        posterior_winner(
          ak     = a,
          dfk    = dfk,
          theta  = theta,
          Mu_dlt = Mu_dlt
        )
      }

      if(count_em == 1){
        init <- theta["mu_a"]
      }else{
        init <- Eak[counter]
      }

      res <- local({
        invisible(capture.output({
          tmp <- MCMC(
            p        = log_post_ak,
            n        = n_mcmc,
            init     = init,
            scale    = 0.1,
            adapt    = TRUE,
            acc.rate = 0.234,
            showProgressBar = FALSE
          )
        }))
        tmp
      })

      Ak[counter,] <- sample(res$samples[-c(1:burnin)],M)

      t <- dfk$t
      dlt_t <- diff(t)
      dlt_Y <- diff(dfk$Y)
      nk <- length(dlt_t)

      mu_mat <- sapply(
        1:nk,
        function(i)
          Mu_dlt(Ak[counter,], ti = t[i+1], ti_1 = t[i])
      )

      res2 <- sweep(mu_mat, 2, dlt_Y, FUN = "-")^2

      E_res2 <- colMeans(res2)  # E[(... )^2]

      num_sigma <- num_sigma + sum(E_res2 / (dlt_t^2))
      den_sigma <- den_sigma + sum(1 / dlt_t)

      counter <- counter + 1
    }

    Eak <- rowMeans(Ak)
    Eak2 <- rowMeans(Ak^2)

    theta["mu_a"] <- mean(Eak)
    theta["sigma2_a"] <- mean(Eak2-2*Eak*theta["mu_a"]+theta["mu_a"]^2)

    theta["sigma2"] <- num_sigma / den_sigma

    dif <- abs(theta - theta_hist[[count_em]])
    error <- max(dif)

    count_em <- count_em + 1
    theta_hist[[count_em]] <- theta

    log_msg(
      sprintf(
        "EM iter %3d | error = %.3e | mu_a = %.4f | sigma2_a = %.4f | sigma2 = %.4f\n",
        count_em,
        error,
        theta["mu_a"],
        theta["sigma2_a"],
        theta["sigma2"]
      )
    )
  }
  return(theta)
}
