#' Generates the density distribution or the accumulated distribution function Exp-G
#'
#' @param G Cumulative distribution function (baseline)
#' @param density By default, `density = FALSE`, that is, it returns the cumulative
#' distribution function \eqn{Exp-G}. If `density = TRUE`, then the \eqn{Exp-G} density
#' function will be reinforced.
#' @return The probability density function or the cumulative distribution
#' function will be returned.
#' @export
#' @details This function generates the probability density function or the
#' accumulated distribution function \eqn{Exp-G}, for any specified \eqn{G}.
#' Exp-G will have one more parameter, being it \eqn{s > 0}.
#' @examples
#' cdf_w <- function(x, beta, lambda) {
#'   1 - exp(-(lambda * x)^beta)
#' }
#'
#' # Exp-Weibull density:
#' pdf_exp_w <- exp_G(G = cdf_w, density = TRUE)
#'
#' # Testing if exp_w is density:
#' integrate(
#'    f = pdf_exp_w,
#'    lower = 0,
#'    upper = Inf,
#'    beta = 1.2,
#'    lambda = 1.3,
#'    s = 2.1
#' )
#'
#' # Exp-Weibull (accumulated distribution):
#' cdf_exp_w <- exp_G(G = cdf_w, density = FALSE)
#'
#' # Testing whether it is cumulative distribution:
#' cdf_exp_w(x = 0, beta = 1.2, lambda = 1.3, s = 2.1) # in 0
#' cdf_exp_w(x = Inf, beta = 1.2, lambda = 1.3, s = 2.1) # in Inf
exp_G <- function(G, density = FALSE) {

  if(density){
    function(x, s, ...) {
      numDeriv::grad(func = function(x) G(x, ...) ^ s, x = x)
    }
  } else {
    function(x, s, ...) {
      G(x = x, ...) ^ s
    }
  }
}

#' Mathematical term used in the sum of equation 19 of the paper
#'
#' @param n Value that indexes a sum of equation 19.
#' @param k Value that indexes a sum of equation 19.
#' @param alpha Additional parameters introduced by the new model proposed in the paper.
#' @param theta Additional parameters introduced by the new model proposed in the paper.
#' @param p_n Zero-truncated power series distribution.
#' @details \eqn{n} and \eqn{k} are indexed values in the sums of equation
#' 19 of the paper.
#' @return Return of expression  \eqn{\upsilon_{n,k}(\alpha, \beta)} ot the paper.
#' @export
upsilon <- function(n, k, alpha, theta, p_n) {
  (-1)^k * choose(-n, k) * alpha^(-n) * (1 - alpha^(-1))^k * p_n(n = n, theta = theta)
}

#' MOTP-Weibull density function (MOTPW)
#'
#' @param x support of the probability density function.
#' @param alpha Parameter.
#' @param theta Parameter.
#' @param beta Parameter.
#' @param lambda Parameter.
#'
#' @examples
#'
#' # Checking the implementation of the theorical density function:
#' integrate(
#'   pdf_theorical,
#'   lower = 0,
#'   upper = 20,
#'   alpha = 1,
#'   theta = 0.5,
#'   beta = 1,
#'   lambda = 1
#' )
#' @export
# Theoretical density (eq. 9 of the paper):
pdf_theorical <- function(x, alpha, theta, beta, lambda) {
   u <- (lambda * x)^beta

   part_1 <- (alpha * theta * beta * lambda^beta * x^(beta - 1) * exp(-u)) /
   ((exp(theta) - 1) * (1 - (1 - alpha) * exp(-u))^2)

   part_2 <- exp(theta * (1 - exp(-u)) / (1 - (1 - alpha) * exp(-u)))

   part_1 * part_2
}

#' Equation 19 of the paper
#' @importFrom future plan
#' @importFrom furrr future_map_dbl
#' @importFrom stats pweibull
#' @param x Support of the probability density function.
#' @param G Baseline cumulative distribution function.
#' @param p_n Zero-truncated power series distribution.
#' @param N Maximum number of sums of the external sum of Eq. 19 of the paper.
#' @param K Maximum number of sums of the interior sum of Eq. 19 of the paper.
#' @param alpha Parameter added by the proposed family.
#' @param theta Parameter added by the proposed family.
#' @param parallel If `parallel = TRUE`, the code will be executed in parallel,
#' taking advantage of the predecessor's multicolors, if any. If `parallel = FALSE`,
#' the code will be executed serially.
#' @param ... Baseline distribution parameters.
#' @export
#' @examples
#' # Implementing the Weibull distribution function (baseline). This function is
#' # used in the eq_19 function. Another baseline G can be passed and tested in
#' # eq_19.
#' cdf_w <- function(x, beta, lambda) {
#'    1 - exp(-(lambda * x)^beta)
#' }
#'
#' pdf_aprox <- function(x, alpha, theta, lambda, beta) {
#'   eq_19(
#'     x = x,
#'     G = cdf_w,
#'     p_n = pf_ztp,
#'     N = 200L,
#'     K = 100L,
#'     alpha = alpha,
#'     theta = theta,
#'     beta = beta,
#'     lambda = lambda
#'   )
#' }
#'
#' # True parameters
#' x = 1.28
#' alpha = 1.2
#' theta = 1.5
#' beta = 1.33
#' lambda = 2
#'
#' # Checking theorical density and approximate density for a finite sum of Eq 19 of the paper
#' pdf_aprox(
#'    x = x,
#'    alpha = alpha,
#'    theta = theta,
#'    lambda = lambda,
#'    beta = beta
#' )
#' pdf_theorical(
#'    x = x,
#'    alpha = alpha,
#'    theta = theta,
#'    beta = beta,
#'    lambda = lambda
#' )
#' @export
#'
eq_19 <- function(x, G, p_n, N = 50L, K = 50L, alpha, theta, parallel = FALSE, ...) {

  pi <- exp_G(G = G, density = TRUE)

  inner_sum <- function(n) {
    step_1 <- function(k) {
      upsilon(
        n = n,
        k = k,
        alpha = alpha,
        theta = theta,
        p_n = p_n
      ) * pi(x = x, s = n + k, ...)
    }
    vapply(
      X = 0L:K,
      FUN = step_1,
      FUN.VALUE = double(length = 1L)
      ) %>% sum()
  }

  if(parallel){
    future::plan("multisession", workers = parallel::detectCores() - 1L)
    furrr::future_map_dbl(.x = 1L:N, .f = inner_sum) %>%
      sum()
  } else {
    vapply(
       X = 1L:N,
       FUN = inner_sum,
       FUN.VALUE = double(length = 1L)
    ) %>% sum()
  }
}
