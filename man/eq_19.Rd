% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exp-G.R
\name{eq_19}
\alias{eq_19}
\title{Equation 19 of the paper}
\usage{
eq_19(x, G, p_n, N = 50L, K = 50L, alpha, theta, parallel = FALSE, ...)
}
\arguments{
\item{x}{Support of the probability density function.}

\item{G}{Baseline cumulative distribution function.}

\item{p_n}{Zero-truncated power series distribution.}

\item{N}{Maximum number of sums of the external sum of Eq. 19 of the paper.}

\item{K}{Maximum number of sums of the interior sum of Eq. 19 of the paper.}

\item{alpha}{Parameter added by the proposed family.}

\item{theta}{Parameter added by the proposed family.}

\item{parallel}{If \code{parallel = TRUE}, the code will be executed in parallel,
taking advantage of the predecessor's multicolors, if any. If \code{parallel = FALSE},
the code will be executed serially.}

\item{...}{Baseline distribution parameters.}
}
\description{
Equation 19 of the paper
}
\examples{
# Implementing the Weibull distribution function (baseline). This function is
# used in the eq_19 function. Another baseline G can be passed and tested in
# eq_19.
cdf_w <- function(x, beta, lambda) {
   1 - exp(-(lambda * x)^beta)
}

pdf_aprox <- function(x, alpha, theta, lambda, beta) {
  eq_19(
    x = x,
    G = cdf_w,
    p_n = pf_ztp,
    N = 200L,
    K = 100L,
    alpha = alpha,
    theta = theta,
    beta = beta,
    lambda = lambda
  )
}

# True parameters
x = 1.28
alpha = 1.2
theta = 1.5
beta = 1.33
lambda = 2

# Checking theorical density and approximate density for a finite sum of Eq 19 of the paper
pdf_aprox(
   x = x,
   alpha = alpha,
   theta = theta,
   lambda = lambda,
   beta = beta
)
pdf_theorical(
   x = x,
   alpha = alpha,
   theta = theta,
   beta = beta,
   lambda = lambda
)
}
