#' Zero-truncated Poisson probability function
#'
#' @param n Function support, with \eqn{n = 1, 2, \cdots}.
#' @param theta Theta parameter, such that \eqn{\theta > 0}.
#'
#' @return Calculated probability (vector)
#' @export
#'
#' @examples
#' sapply(
#'    X = 1L:1000L,
#'    FUN = function(i)
#'    pf_ztp(n = i, theta = 1.2)
#' ) %>% sum
pf_ztp <- function(n, theta) {
  theta ^ n / ((exp(theta) - 1) * factorial(n))
}

#' Logarithmic probability function
#' @param n Function support, with \eqn{n = 1, 2, \cdots}.
#' @param theta A probability value, i.e, \eqn{0 < \theta < 1}.
#'
#' @return Calculated probability (vector)
#' @export
#'
#' @examples
#' sapply(
#'      X = 1L:1000L,
#'      FUN = function(i)
#'      pf_logarithmic(n = i, theta = 0.5)
#' ) %>% sum
pf_logarithmic <- function(n, theta) {

  f <- function(theta) {
    -theta^n / (n * log(1 - theta))
  }

  if(n == 0) {
    (1 - theta) + theta * f(0)
  } else {
    f(theta)
  }
}

#' Geometric probability function
#' @param n Function support, with \eqn{n = 1, 2, \cdots}.
#' @param theta A probability value, i.e, \eqn{0 < \theta < 1}.
#'
#' @return Calculated probability (vector)
#' @export
#'
#' @examples
#' sapply(
#'    X = 1L:1000L,
#'    FUN = function(i)
#'    pf_geometric(n = i, theta = 0.5)
#' ) %>% sum
pf_geometric <- function(n, theta) {
   theta * (1 - theta)^(n - 1)
}
