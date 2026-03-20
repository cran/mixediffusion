\name{datasim01}
\alias{datasim01}
\docType{data}
\title{
Simulated diffusion dataset 01
}
\description{
A simulated dataset consisting of discretely observed trajectories generated
from a one-dimensional Wiener diffusion model with random effects in the drift
term.
}
\usage{
data(datasim01)
}
\format{
A data frame with observations from multiple experimental units and the
following variables:
\describe{
  \item{\code{t}}{Observation time.}
  \item{\code{unit}}{Unit identifier.}
  \item{\code{Y}}{Observed process value.}
}
}
\details{
The data were generated from the stochastic differential equation
\deqn{
dY_k(t) = a_k \, t \, dt + \sigma \, dW_k(t),
}
where the unit-specific random effects follow a normal distribution
\deqn{
a_k \sim \mathcal{N}(10, 3),
}
and the diffusion variance is given by \eqn{\sigma^2 = 500}.

The process was simulated for \eqn{K = 50} units, each observed at
\eqn{n = 200} equally spaced time points, with initial condition
\eqn{Y_{k0} = 1}.
}
\source{
Simulated data generated for illustrative purposes.
}
\value{
A data frame with simulated data.
}
\references{
None.
}
\examples{
data(datasim01)
plot_paths(df = datasim01)
}
\keyword{datasets}
