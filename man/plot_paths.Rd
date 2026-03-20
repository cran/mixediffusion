\name{plot_paths}
\alias{plot_paths}
\title{
Plot panel trajectories
}
\description{
Plots multiple trajectories from panel or longitudinal data grouped by unit.
}
\usage{
plot_paths(df,
  col = NULL,
  lwd = 1,
  xlab = "t",
  ylab = "Y",
  main = NULL,
  ...)
}
\arguments{
  \item{df}{
  Data frame containing the observed data. It must include the columns
  \code{t} (time), \code{unit} (unit identifier) and \code{Y} (observed values).
  }

  \item{col}{
  Optional vector of colors, one per unit. If \code{NULL}, a grayscale palette
  is used.
  }

  \item{lwd}{
  Line width used for the trajectories.
  }

  \item{xlab}{
  Label for the x-axis.
  }

  \item{ylab}{
  Label for the y-axis.
  }

  \item{main}{
  Optional main title for the plot.
  }

  \item{...}{
  Additional graphical parameters passed to \code{plot}.
  }
}

\value{
A ggplot graph.
}
\examples{
data(datasim02)
plot_paths(datasim02)
}
