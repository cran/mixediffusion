make_Mu_cond <- function(media) {

  expr <- gsub("\\s+", "", media)

  Mu <- switch(
    expr,
    "a"       = Mu_cond_const,
    "at^1"    = Mu_cond_at1,
    "at^2"    = Mu_cond_at2,
    "exp(at)" = Mu_cond_exp,
    "cos(at)" = Mu_cond_cos,
    "sin(at)" = Mu_cond_sin,
    NULL
  )

  if (is.null(Mu)) {
    stop(
      sprintf(
        "Drift '%s' not recognized. Valid options are: a, at^1, at^2, exp(at), cos(at), sin(at).",
        media
      )
    )
  }

  return(Mu)
}
