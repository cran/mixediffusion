\name{datasim02}
\alias{datasim02}
\docType{data}
\title{
Simulated diffusion dataset 02
}
\description{
A simulated dataset consisting of discretely observed trajectories generated
from a one-dimensional Wiener diffusion model with random effects in the drift
term.
}
\usage{
data(datasim02)
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
dY_k(t) = a_k \sin(\pi t)\,dt + \sigma\,dW_k(t),
}
where the unit-specific random effects follow a normal distribution
\deqn{
a_k \sim \mathcal{N}(10, 3),
}
and the diffusion variance is given by \eqn{\sigma^2 = 50}.

The integrated drift used for simulation is given by
\deqn{
\mu_{\Delta}(a_k, t_i, t_{i-1})
=
\frac{a_k}{\pi}
\left\{
\cos(\pi t_{i-1}) - \cos(\pi t_i)
\right\}.
}

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
data(datasim02)
plot_paths(df = datasim02)
}
\keyword{datasets}
