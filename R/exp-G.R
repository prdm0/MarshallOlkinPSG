#' Generates the density distribution or the accumulated distribution function Exp-G
#'
#' @param G Cumulative distribution function (baseline)
#' @param density By default, `density = FALSE`, that is, it returns the cumulative
#' distribution function \eqn{Exp-G}. If `density = TRUE`, then the \eqn{Exp-G} density
#' function will be reinforced.
#'
#' @return The probability density function or the cumulative distribution
#' function will be returned.
#' @export
#' @details  This function generates the probability density function or the
#' accumulated distribution function \eqn{Exp-G}, for any specified \eqn{G}.
#' Exp-G will have one more parameter, being it \eqn{s > 0}.
#' @examples
#' cdf_w <- function(x, alpha, beta) {
#'    pweibull(q = x, shape = alpha, scale = beta)
#' }
#'
#' # Exp-Weibull density:
#' pdf_exp_w <- exp_G(G = pdf_w, density = TRUE)
#'
#' # Testing if exp_w is density:
#' integrate(
#'    f = pdf_exp_w,
#'    lower = 0,
#'    upper = Inf,
#'    alpha = 1.2,
#'    beta = 1.3,
#'    s = 2.1
#' )
#'
#' # Exp-Weibull (accumulated distribution):
#' cdf_exp_w <- exp_G(G = pdf_w, density = FALSE)
#'
#' # Testing whether it is cumulative distribution:
#' cdf_exp_w(x = 0, alpha = 1.2, beta = 1.3, s = 2.1) # in 0
#' cdf_exp_w(x = Inf, alpha = 1.2, beta = 1.3, s = 2.1) # in Inf
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
#' @details \eqn{n} and \eqn{k} are indexed values in the sums of equation
#' 19 of the paper.
#' @return Return of expression  \eqn{\upsilon_{n,k}(\alpha, \beta)} ot the paper.
#' @export
upsilon <- function(n, k, alpha, theta, p_n) {
  (-1)^k * choose(-n, k) * alpha^(-n) * (1 - alpha^(-1))^k * p_n(n = n, theta = theta)
}

#' Equation 19 of the paper
#'
#' @examples
#' cdf_w <- function(x, a, b) {
#'    pweibull(q = x, shape = a, scale = b)
#' }
#' eq_19(alpha, theta, p_n, G)
#' @export
eq_19 <- function(x, G, p_n, N = 100L, K = 100L, alpha, theta, ...) {

  pi <- exp_G(G = G, density = TRUE)

  inner_sum <- function(n) {
    step_1 <- function(k) {
      upsilon(
        n = n,
        k = k,
        alpha = alpha,
        theta = theta,
        p_n = p_n
      ) * pi(s = n + k, ...)
    }
    vapply(
      X = 0L:K,
      FUN = step_1,
      FUN.VALUE = double(length = 1L)
    ) %>% sum()
  }

  vapply(
    X = 1L:N,
    FUN = inner_sum,
    FUN.VALUE = double(length = 1L)
  ) %>% sum()
}

