library(dplyr)
library(tidyr)
library(ggplot2)

tbl_duration = readRDS("duration_samples.rds")
logit_hazard_mean = readRDS("logit_hazard_mean.rds")
logit_hazard_cov = readRDS("logit_hazard_cov.rds")

logit = function(x) log(x) - log(1 - x)
expit = function(x) 1 / (1 + exp(-x))

# Happens a negligiable amount!
tbl_duration |>
    filter(time <= 0, f > 0 | F > 0 | S < 1)

tbl_duration |>
    group_by(time) |>
    tidybayes::median_qi(F) |>
    print(n=50)

tbl_duration |>
    tidyr::pivot_longer(!c(.draw, time)) |>
    group_by(time, name) |>
    tidybayes::median_qi(.width = c(0.5, 0.8, 0.95)) |>
    ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
    tidybayes::geom_lineribbon() +
    facet_wrap(~name) +
    scale_x_continuous(breaks = (0:100) * 14, minor_breaks = (0:100) * 2)

logit_hazard_matrix = tbl_duration |>
    filter(between(time, 0, 40)) |>
    mutate(lambda = if_else(lambda <= 0, (1e-7)/2, lambda)) |>
    assertr::verify(lambda > 0 & lambda < 1) |>
    pivot_wider(id_cols = .draw, values_from = lambda, names_from = time) |>
    select(!.draw) |>
    as.matrix() |>
    logit()
stopifnot(all(is.finite(logit_hazard_matrix)))

logit_hazard_cov = logit_hazard_matrix |>
    cov()
stopifnot(all(is.finite(logit_hazard_cov)))

logit_hazard_mean = logit_hazard_matrix |>
    colMeans()


logit_hazard_matrix2 |>
    cov() |>
    diag() |>
    cbind(diag(logit_hazard_cov))
stopifnot(all(is.finite(logit_hazard_cov)))

logit_hazard_matrix2 |>
    colMeans()

# Means from the two methods similar
# Variance quite a bit higher using beta method which is probability good: otherwise unrealistically low!
# Only matters for first 3 elements, and we discard one of these

# Check correlations
logit_hazard_matrix2 |>
    cor() |>
    magrittr::extract(1:6, 1:6)
logit_hazard_cov |>
    cov2cor() |>
    magrittr::extract(1:6, 1:6)
# 0 and 1 quite a bit more correlated on the latter (with each other)
# Generally little difference though

## What do marginals for day 1 (element 2) look like?
# Orig: very, very small (all <= 1e-6)
qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), logit_hazard_mean[2], sqrt(logit_hazard_cov[2,2])) |> expit()
# Updated, slightly bigger but not much (largest 7e-6)
qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), mean(logit_hazard_matrix2[,2]), sd(logit_hazard_matrix2[,2])) |> expit()

## What do marginals for day 2 (element 3) look like?
# Orig: very, very small (all <= 1e-6)
qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), logit_hazard_mean[3], sqrt(logit_hazard_cov[3,3])) |> expit()
# Updated, slightly bigger but not much (largest 1e-4)
qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), mean(logit_hazard_matrix2[,3]), sd(logit_hazard_matrix2[,3])) |> expit()