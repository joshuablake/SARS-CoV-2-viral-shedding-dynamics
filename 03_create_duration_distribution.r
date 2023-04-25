library(dplyr)
library(rstan)
library(tidyr)

n_param_samples = 1e3
n_indiv_per_param = 1e6
n_group = 1 # Group number to calculate distribution for (see 01_run_RNA.R)

logit = function(x) log(x) - log(1 - x)

logsumexp <- function(a, b) {
    max_x <- pmax(a, b)  # find the maximum value in x
    sum_exp <- exp(a - max_x) + exp(b - max_x)  # subtract max_x to avoid overflow
    log_sum_exp <- log(sum_exp) + max_x  # add max_x back in to get the log-sum-exp
    return(log_sum_exp)
}

calc_durations = function(theta) {
    term1 = (exp(theta[2,]) + exp(theta[3,])) * (theta[1,] + logsumexp(theta[2,], theta[3])) / exp(theta[2,] + theta[3,])
    term2 = theta[2,] * exp(-theta[3,])
    term3 = theta[3,] * exp(-theta[2,])
    return(term1 - term2 - term3)
}

fit = readRDS(here::here("fit2.rds"))
max_iterations = floor((fit@sim$iter - fit@sim$warmup) * fit@sim$chains / fit@sim$thin)
draw_nums = sample.int(max_iterations, n_param_samples)
draws = rstan::extract(fit, c("vgrow", "v_sd", "Lc"))

# Create distribution for each parameter
tbl_duration = purrr::map_dfr(
    draw_nums,
    function(.draw) {
        # Simulate three standard normals per individual (one per random effect)
        z = rnorm(n_indiv_per_param * 3) |>
            matrix(nrow = 3)

        # Parameters determining conversion
        mu = draws$vgrow[.draw, , n_group]
        C = draws$Lc[.draw,,]
        delta = draws$v_sd[.draw,]

        # Convert simulated normals to correct distribution
        theta = mu + delta * (C %*% z)
        # For each parameter set, calculate the duration
        duration = calc_durations(theta)
        stopifnot(all(duration >= 0))
        stopifnot(length(duration) == n_indiv_per_param)
        return(
            tibble::tibble(
                .draw = .draw,
                t = 0:50,
                f = purrr::map_dbl(t, ~mean(duration <= .x + 0.5 & duration > .x - 0.5)),
                F = purrr::map_dbl(t, ~mean(duration <= .x + 0.5)),
            )
        )
    },
    .progress = TRUE
) |>
    mutate(
        S = c(1, 1 - lag(F)[-1]),
        lambda = if_else(S == 0, 0, f / S),
    )

tbl_duration |>
    group_by(t) |>
    tidybayes::median_qi(F) |>
    print(n=50)

tbl_duration |>
    tidyr::pivot_longer(!c(.draw, t)) |>
    group_by(t, name) |>
    tidybayes::median_qi(.width = c(0.5, 0.8, 0.95)) |>
    ggplot(aes(t, value, ymin = .lower, ymax = .upper)) +
    tidybayes::geom_lineribbon() +
    facet_wrap(~name) +
    scale_x_continuous(breaks = (0:100) * 14, minor_breaks = (0:100) * 2)

logit_hazard_matrix = tbl_duration |>
    filter(between(t, 0, 40)) |>
    mutate(lambda = if_else(lambda <= 0, (1e-6)/2, lambda)) |>
    assertr::verify(lambda > 0 & lambda < 1) |>
    pivot_wider(id_cols = .draw, values_from = lambda, names_from = t) |>
    select(!.draw) |>
    as.matrix() |>
    logit()
stopifnot(all(is.finite(logit_hazard_matrix)))

logit_hazard_cov = logit_hazard_matrix |>
    cov()
stopifnot(all(is.finite(logit_hazard_cov)))

logit_hazard_mean = logit_hazard_matrix |>
    colMeans()

saveRDS(tbl_duration, "duration_samples.rds")
saveRDS(logit_hazard_mean, "logit_hazard_mean.rds")
saveRDS(logit_hazard_cov, "logit_hazard_cov.rds")
