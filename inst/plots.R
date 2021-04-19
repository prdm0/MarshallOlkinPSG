# remotes::install_github("wilkelab/cowplot")

library(cowplot)
library(ggplot2)
library(magrittr)

p1 <- eval_plot_moptw(
  n_points = 200L,
  N = 1L,
  K = 1L,
  xmin = 0.01,
  xmax = 2,
  alpha = 1.2,
  theta = 1.5,
  beta = 1.33,
  lambda = 2
)

p2 <- eval_plot_moptw(
  n_points = 200L,
  N = 2L,
  K = 2L,
  xmin = 0.01,
  xmax = 2,
  alpha = 1.2,
  theta = 1.5,
  beta = 1.33,
  lambda = 2
)

p3 <- eval_plot_moptw(
  n_points = 200L,
  N = 3L,
  K = 3L,
  xmin = 0.01,
  xmax = 2,
  alpha = 1.2,
  theta = 1.5,
  beta = 1.33,
  lambda = 2
)

p4 <- eval_plot_moptw(
  n_points = 200L,
  N = 4L,
  K = 4L,
  xmin = 0.01,
  xmax = 2,
  alpha = 1.2,
  theta = 1.5,
  beta = 1.33,
  lambda = 2
)

plot_grid(p1, p2, p3, p4, labels = LETTERS[1L:4L], label_size = 12)
  ggsave("~/Downloads/MarshallOlkinPSG/inst/numerical_evaluation.png", width = 8.4, height = 8.4)
