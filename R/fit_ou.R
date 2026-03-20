fit_ou <- function(df,mu = "at^1",tol = 1e-4,
                  max_iter = 100,theta = NULL,
                  M = 100,verbose = TRUE,
                  mu_cond = NULL,n_mcmc = 1000, burnin = 500){
  if(!is.null(mu_cond)){
    Mu_cond <- mu_cond
  }else{
    Mu_cond <- make_Mu_cond(mu)
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
      sigma2   = 1,
      lambda = 0.1
    )
  }
  theta_hist <- list()
  theta_hist[[1]] <- theta
  K <- length(unique(df$unit))
  #M <- 100
  #tol <- 1e-4
  count_em <- 1
  error <- 1e300
  #max_iter <- 100
  while (error > tol & count_em < max_iter){
    Ak <- matrix(ncol = M, nrow = K)
    counter <- 1
    for(k in unique(df$unit)){
      dfk <- df %>%
        filter(unit == k)
      log_post_ak <- function(a){
        posterior_ou(
          ak     = a,
          dfk    = dfk,
          theta  = theta,
          Mu_cond = Mu_cond
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
      counter <- counter + 1
    }
    Eak <- rowMeans(Ak)
    Eak2 <- rowMeans(Ak^2)
    theta["mu_a"] <- mean(Eak)
    theta["sigma2_a"] <- mean(Eak2-2*Eak*theta["mu_a"]+theta["mu_a"]^2)

    opt <- optim(
      par     = theta["lambda"],
      fn      = Q_lambda,
      method  = "Brent",
      lower   = 1e-7,
      upper   = 50,
      control = list(fnscale = -1),
      df = df,
      Ak_draws = Ak,
      sigma2 = theta["sigma2"],
      Mu_cond = Mu_cond
    )

    theta["lambda"] <- opt$par

    theta["sigma2"] <- update_sigma2(lambda = theta["lambda"],
                                     df = df,
                                     Ak_draws = Ak,
                                     Mu_cond = Mu_cond)


    dif <- abs(theta - theta_hist[[count_em]])
    error <- max(dif)
    #print(error)
    count_em <- count_em + 1
    theta_hist[[count_em]] <- theta
    log_msg(
      sprintf(
        "EM iter %3d | error = %.3e | mu_a = %.4f | sigma2_a = %.4f | sigma2 = %.4f | lambda = %.4f\n",
        count_em,
        error,
        theta["mu_a"],
        theta["sigma2_a"],
        theta["sigma2"],
        theta["lambda"]
      )
    )
  }
  return(theta)
}
