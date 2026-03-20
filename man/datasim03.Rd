\name{datasim03}
\alias{datasim03}
\docType{data}
\title{
Simulated diffusion dataset 03
}
\description{
A simulated dataset consisting of discretely observed trajectories generated
from a one-dimensional Ornstein--Uhlenbeck diffusion model with a time-dependent
mean function and random effects.
}
\usage{
data(datasim03)
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
dY_k(t) = \lambda \{ a_k t - Y_k(t) \} \, dt + \sigma \, dW_k(t),
}
where the unit-specific random effects follow a normal distribution
\deqn{
a_k \sim \mathcal{N}(10, 3),
}
the mean-reversion parameter is given by \eqn{\lambda = 0.75}, and the diffusion
variance satisfies \eqn{\sigma^2 = 500}.

The process was simulated for \eqn{K = 30} units, each observed at
\eqn{n = 200} equally spaced time points over the interval \eqn{[0, 10]},
with initial condition \eqn{Y_k(0) = 0}.
}
\source{
Simulated data generated for illustrative purposes.
}
\references{
None.
}
\value{
A data frame with simulated data.
}
\examples{
data(datasim03)
plot_paths(df = datasim03)
}
\keyword{datasets}
