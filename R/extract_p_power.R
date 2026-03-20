extract_p_power <- function(expr) {

  expr <- gsub("\\s+", "", expr)

  # caso t^(...)
  m1 <- regexpr("t\\^\\([^\\)]+\\)", expr)

  if (m1 != -1) {
    p_str <- regmatches(expr, m1)
    p_str <- sub("t\\^\\(", "", p_str)
    p_str <- sub("\\)", "", p_str)

    p <- tryCatch(eval(parse(text = p_str)),
                  error = function(e) NA)

    if (is.na(p)) {
      stop("No se pudo evaluar la potencia p")
    }

    return(as.numeric(p))
  }

  # caso t^p directo
  m2 <- regexpr("t\\^[-+]?[0-9]*\\.?[0-9]+", expr)

  if (m2 == -1) {
    stop("No se pudo detectar una potencia del tipo t^p")
  }

  p <- as.numeric(sub("t\\^", "", regmatches(expr, m2)))

  return(p)
}
