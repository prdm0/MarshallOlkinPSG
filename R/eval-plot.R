# # Theoretical density (eq. 9 of the paper):
# pdf_motpw <- function(x, alpha, theta, beta, lambda) {
#   u <- (lambda * x)^beta
#
#   part_1 <- (alpha * theta * beta * lambda^beta * x^(beta - 1) * exp(-u)) /
#     ((exp(theta) - 1) * (1 - (1 - alpha) * exp(-u))^2)
#
#   part_2 <- exp(theta * (1 - exp(-u)) / (1 - (1 - alpha) * exp(-u)))
#
#   part_1 * part_2
# }
#
# pdf_motpw(x = 0.1, alpha = 0.5, theta = 0.5, beta = 0.5, lambda = 0.5)
#
# cdf_w <- function(x, beta, lambda) {
#   pweibull(q = x, shape = beta, scale = lambda)
# }
#
# pdf_motpw_aprox <- function(x) {
#   eq_19(
#     x = x,
#     G = cdf_w,
#     p_n = pf_ztp,
#     N = 200L,
#     K = 200L,
#     alpha = 1.5,
#     theta = 1.5,
#     beta = 1.5,
#     lambda = 1.5
#   )
# }

#' Plot Eval MOPTW
#' @export
#' @examples
#' eval_plot_moptw(
#'    n_points = 5L,
#'    N = 100L,
#'    K = 100L,
#'    xmin = 0.01,
#'    xmax = 2,
#'    alpha = 0.5,
#'    theta = 0.5,
#'    beta = 0.5,
#'    lambda = 0.5
#' )
eval_plot_moptw <- function(n_points = 5L, N = 50L, K = 50L, xmin, xmax, alpha, theta, beta, lambda) {
  seq_x <- seq(xmin, xmax, length.out = n_points)

  cdf_w <- function(x, beta, lambda) {
    pweibull(q = x, shape = beta, scale = lambda)
  }

  # Theoretical density (eq. 9 of the paper):
  pdf_motpw <- function(x, alpha, theta, beta, lambda) {
    u <- (lambda * x)^beta

    part_1 <- (alpha * theta * beta * lambda^beta * x^(beta - 1) * exp(-u)) /
      ((exp(theta) - 1) * (1 - (1 - alpha) * exp(-u))^2)

    part_2 <- exp(theta * (1 - exp(-u)) / (1 - (1 - alpha) * exp(-u)))

    part_1 * part_2
  }
  y_theorical <-
    pdf_motpw(
       x = seq_x,
       alpha = alpha,
       theta = theta,
       beta = beta,
       lambda = lambda
    )

  y_aprox_step <- function(x) {
    eq_19(
       x = x,
       G = cdf_w,
       p_n = pf_ztp,
       N = N,
       K = K,
       alpha = alpha,
       theta = theta,
       beta = beta,
       lambda = lambda
    )
  }

  y_aprox <- vapply(X = seq_x, FUN = y_aprox_step, FUN.VALUE = double(length = 1L))

  df <- data.frame(x = rep(seq_x, 2L), y = c(y_theorical, y_aprox),
        group = c(rep("theorical", length(y_theorical)), rep("aprox", length(y_aprox))))

  ggplot2::ggplot(df, ggplot2::aes(x, y, group = group)) +
     ggplot2::geom_line()

}
