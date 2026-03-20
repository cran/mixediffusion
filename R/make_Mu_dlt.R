make_Mu_dlt <- function(media) {

  expr <- gsub("\\s+", "", media)

  # --- familia potencia
  if (grepl("^at\\^", expr)) {

    p <- extract_p_power(expr)

    return(
      function(ak, ti, ti_1) {
        Mu_dlt_power(ak, ti, ti_1, p)
      }
    )
  }

  # --- familias cerradas exactas
  Mu <- switch(
    expr,
    "exp(at)" = Mu_dlt_exp,
    "cos(at)" = Mu_dlt_cos,
    "sin(at)" = Mu_dlt_sin,
    NULL
  )

  if (is.null(Mu)) {
    stop(
      sprintf(
        "Drift '%s' not recognized. Valid options are: exp(at), sin(at), at^p",
        media
      )
    )
  }

  return(Mu)
}
