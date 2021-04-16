#' Plot Eval MOPTW
#' @importFrom ggplot2 ggplot geom_line aes
#' @param n_points Number of points to form the plot;
#' @param N Maximum number of sums of the external sum of Eq. 19 of the paper.
#' @param K Maximum number of sums of the interior sum of Eq. 19 of the paper.
#' @param xmin Minimum value for the domain. This parameter allows you to control the x axis
#' @param xmax Maximum value for the domain. This parameter allows you to control the x axis
#' @param parallel If `parallel = TRUE`, the code will be executed in parallel, taking advantage of the predecessors multicolors, if any. If `parallel = FALSE`, the code will be executed serially.
#' @param alpha Parameter added by the proposed family.
#' @param theta Parameter added by the proposed family.
#' @param beta Scale parameter of the Weibull distribution.
#' @param lambda Shape parameter of the Weibull distribution.
#'
#' @return Line graph comparing the theoretical distribution of Eq. 9 of the paper
#' with an approximation, considering finite sums, by Eq. 19 of the paper.
#'
#' @export
#'
#' @examples
#' eval_plot_moptw(
#'    n_points = 2L,
#'    N = 200L,
#'    K = 100L,
#'    xmin = 0.01,
#'    xmax = 2,
#'    alpha = 1.2,
#'    theta = 1.5,
#'    beta = 1.33,
#'    lambda = 2
#' )
eval_plot_moptw <- function(n_points = 5L, N = 200L, K = 100L, xmin, xmax, alpha,
                            theta, beta, lambda, parallel = FALSE) {
  seq_x <- seq(xmin, xmax, length.out = n_points)

  cdf_w <- function(x, beta, lambda) {
    1 - exp(-(lambda * x)^beta)
  }

  # Theoretical density (eq. 9 of the paper):
  pdf_motpw <- function(x, alpha, theta, beta, lambda) {
    u <- (lambda * x)^beta

    part_1 <- (alpha * theta * beta * lambda^beta * x^(beta - 1) * exp(-u)) /
      ((exp(theta) - 1) * (1 - (1 - alpha) * exp(-u))^2)

    part_2 <- exp(theta * (1 - exp(-u)) / (1 - (1 - alpha) * exp(-u)))

    part_1 * part_2
  }

  y_aprox_step <- function(x) {
    eq_19(
       x = x,
       G = cdf_w,
       p_n = pf_ztp,
       N = N,
       K = K,
       alpha = alpha,
       theta = theta,
       lambda = lambda,
       beta = beta,
       parallel = parallel
    )
  }

  y_theorical <-
    pdf_motpw(
      x = seq_x,
      alpha = alpha,
      theta = theta,
      lambda = lambda,
      beta = beta
    )

  y_aprox <- vapply(X = seq_x, FUN = y_aprox_step, FUN.VALUE = double(length = 1L))

  df <- data.frame(x = rep(seq_x, 2L), y = c(y_theorical, y_aprox),
        name = c(rep("theorical", length(y_aprox)), rep("aprox", length(y_theorical))))

  df %>%
     ggplot(aes(x = x, y = y, color = name)) +
     geom_line()

  df

}
