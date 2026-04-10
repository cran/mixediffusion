\name{fit_ou}
\alias{fit_ou}
\title{
Inference about the parameters of an Ornstein--Uhlenbeck process
}
\description{
Implements the EM algorithm to perform inference on the parameters of an
Ornstein--Uhlenbeck process with mixed drift effects.}
\usage{
fit_ou(df,
  mu = "at^1",
  tol = 1e-4,
  max_iter = 100,
  theta = NULL,
  M = 100,
  verbose = TRUE,
  mu_cond = NULL,
  n_mcmc = 1000,
  burnin = 500)
}
\arguments{
  \item{df}{Data frame with the observed data. It must include the columns
  \code{t} (time), \code{unit} (unit identifier) and \code{Y} (the observed
  trajectory).}

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

  \item{mu_cond}{Optional user-supplied function defining the conditional mean of
  the process. If \code{NULL}, it is constructed automatically from \code{mu}.}
  \item{n_mcmc}{Number of MCMC iterations used in the E-step to sample the random effects.}
  \item{burnin}{Number of initial MCMC iterations discarded.}
}
\details{
The model is a one-dimensional Ornstein--Uhlenbeck diffusion defined by
\deqn{
dY_{kt} = \lambda\{\mu(t, a_k) - Y_{kt}\}\,dt + \sigma\,dW_{kt},
}

where \eqn{\mu(t, a_k)} is a user-specified drift function depending on a
unit-specific random effect \eqn{a_k \sim \mathcal{N}(\mu_a, \sigma_a^2)}.


The parameter \eqn{\lambda} controls the strength of mean reversion toward the
time-dependent mean.

For discretely observed trajectories, the conditional mean is given by
\deqn{
E\{Y_{kt_i} \mid Y_{kt_{i-1}}, a_k\}
=
Y_{kt_{i-1}} e^{-\lambda \Delta t_i}
+
\lambda \int_{t_{i-1}}^{t_i}
e^{-\lambda (t_i - s)} \mu(s,a_k)\,ds.
}

The function \code{mu_cond} represents this conditional mean. If not provided, it
is constructed automatically from the drift specification supplied through
\code{mu} using closed-form expressions.
}
\value{
A named numeric vector containing the estimated model parameters:
\item{mu_a}{Estimated mean of the random effects distribution.}
\item{sigma2_a}{Estimated variance of the random effects distribution.}
\item{sigma2}{Estimated diffusion variance of the process.}
\item{lambda}{Estimated mean reversion parameter.}
}


\examples{
library(mixediffusion)
data(datasim03)
plot_paths(df = datasim03)
fit <- fit_ou(df = datasim03, mu = "at^1",
              verbose = FALSE, max_iter = 1)
fit
}


