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
#' @importFrom actuar dlogarithmic
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
  actuar::dlogarithmic(x = n, prob = theta, log = FALSE)
}

#' Geometric probability function
#' @importFrom stats dgeom
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
  dgeom(x = n, prob = theta)
}
