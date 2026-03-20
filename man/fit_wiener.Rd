\name{fit_wiener}
\alias{fit_wiener}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference about the parameters of a Wiener process
}
\description{
Implement the EM algorithm to perform inference on the parameters of a Wiener process with mixed drift effects.}
\usage{
fit_wiener(df,
  mu = "at^1",
  tol = 1e-4,
  max_iter = 100,
  theta = NULL,
  M = 100,
  verbose = TRUE,
  mu_dlt = NULL,
  n_mcmc = 1000,
  burnin = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{Data frame with the observed data. It must include the columns \code{t} (time), \code{unit} (unit identifier) and \code{Y} (the observed trajectory).}
  \item{mu}{Functional form of the drift. Supported drifts include
  \code{"at^p"}, \code{"exp(at)"}, \code{"cos(at)"} and \code{"sin(at)"}.
    The default is \code{"at^1"}.}
  \item{tol}{Convergence tolerance for the EM algorithm. The algorithm stops when
    the maximum absolute difference between the parameter estimates at two
    consecutive EM iterations is smaller than \code{tol}.}
  \item{max_iter}{Maximum number of EM iterations.}
  \item{theta}{Optional named vector of initial parameter values. If \code{NULL},
    default initial values are used.}
  \item{M}{Number of Monte Carlo samples used to approximate the conditional
    expectations in the E-step.}
  \item{verbose}{Logical indicating whether to print EM iteration progress.}

  \item{mu_dlt}{Optional user-supplied function defining the integrated drift
    term. If \code{NULL}, it is constructed automatically from \code{mu}.}
  \item{n_mcmc}{Number of MCMC iterations used in the E-step to sample the random effects.}
  \item{burnin}{Number of initial MCMC iterations discarded.}

}
\details{
The model is a one-dimensional Wiener diffusion defined by
\deqn{
dY_{kt} = \mu(t, a_k)dt + \sigma dW_{kt}},

where \eqn{\mu(t, a_k)} is a user-specified drift function depending on a
unit-specific random effect \eqn{a_k \sim \mathcal{N}(\mu_a, \sigma_a^2)}.

For discretely observed trajectories, the mean of the increments is given by the
integrated drift
\deqn{
\mu_{\Delta}(a_k, t_i, t_{i-1}) = \int_{t_{i-1}}^{t_i} \mu(s, a_k)\, ds.
}

The function \code{mu_dlt} represents this integrated drift. If not provided, it
is constructed automatically from \code{mu} using closed-form expressions for
common drift specifications.
}

\value{
A named numeric vector containing the estimated model parameters:
\item{mu_a}{Estimated mean of the random effects distribution.}
\item{sigma2_a}{Estimated variance of the random effects distribution.}
\item{sigma2}{Estimated diffusion variance of the Wiener process.}
}



\examples{
library(mixediffusion)
data(datasim01)
plot_paths(df = datasim01)
\dontrun{
fit <- fit_wiener(df = datasim01, mu = "at^1",
                  verbose = FALSE, max_iter = 20)
fit
}

# mu(ak,t) = ak*sin(pi*t)
data(datasim02)
plot_paths(df = datasim02)
mu_dlt_new <- function(ak,ti,ti_1){
  value <- -ak*(cos(pi*ti) - cos(pi*ti_1))
  return(value)
}

fit <- fit_wiener(df = datasim02, mu_dlt = mu_dlt_new,
                  verbose = FALSE, max_iter = 2)
fit
}
