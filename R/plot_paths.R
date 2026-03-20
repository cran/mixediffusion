plot_paths <- function(df,
                          col = NULL,
                          lwd = 1,
                          xlab = "t",
                          ylab = "Y",
                          main = NULL,
                          ...) {

  stopifnot(all(c("t", "Y", "unit") %in% names(df)))

  units <- unique(df$unit)

  if (is.null(col)) {
    col <- gray.colors(length(units), start = 0.3, end = 0.8)
  }

  plot(
    NA,
    xlim = range(df$t),
    ylim = range(df$Y),
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  for (i in seq_along(units)) {
    idx <- df$unit == units[i]
    lines(df$t[idx], df$Y[idx], col = col[i], lwd = lwd)
  }

  invisible(NULL)
}
