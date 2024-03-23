library(dplyr)
library(rstan)

set.seed(123)

data <- readRDS(file="Data/trajectories.RDS")  |>
    as_tibble() |>
    mutate(group = case_when(
        (vaccinated == FALSE) & (WGS == "Pre-Alpha") ~ 1,
        (vaccinated == FALSE) & (WGS == "Alpha") ~ 1,
        (vaccinated == FALSE) & (WGS == "Delta") ~ 1,
        (vaccinated == TRUE) & (WGS == "Delta") ~ 2
    ))  %>%
    filter(!is.na(group), group == 1, !is.na(copy)) |>
    mutate(
        positive = copy > 1,
        log_vl = log(copy) - 3.43,
        indiv_id = dense_rank(id_sub),
    ) |>
    assertr::verify(log_vl > 0 | !positive)

# Copy number
options(mc.cores = 8)
rstan_options(auto_write = TRUE)

init_pos = function() {
    max_measurements = data |>
        summarise(
            t = day[which.max(log_vl)],
            vl = max(log_vl),
            .by = indiv_id
        ) |>
        arrange(indiv_id)
    first_pos = data |>
        filter(positive) |>
        summarise(
            .by = indiv_id,
            t = min(day),
            vl = log_vl[which.min(day)],
        )
    tibble::lst(
        peak_log_viral_load_mean = mean(max_measurements$vl) + rnorm(1),
        peak_log_viral_load_sd = rexp(1, 1/sd(max_measurements$vl)),
        peak_log_viral_load_raw = max_measurements$vl / sd(max_measurements$vl) - mean(max_measurements$vl),
        absolute_t_peak = pmax(max_measurements$t + rnorm(nrow(max_measurements)), first_pos$t + 1),
        logt_peak = log(absolute_t_peak - first_pos$t + first_pos$vl / 5) + rnorm(nrow(first_pos)),
        logt_peak_mean = mean(logt_peak) + rnorm(1),
        logt_peak_sd = rexp(1, 1/sd(logt_peak)),
    )
}

# Pass data to model
data.stan<-list(
    lod = 3.43,
    max_day = 20,
    n_test = nrow(data),
    n_individuals = max(data$indiv_id),
    individual = data$indiv_id,
    day = data$day,
    is_pos = data$positive,
    log_viral_load = data$log_vl
)


fit = stan(
    "original_analysis/adjust_for_Hakki.stan",
    chains=8,
    cores=8,
    data=data.stan,
    control = list(adapt_delta = 0.9999, max_treedepth = 20),
    init = init_pos
)

saveRDS(fit, file = "fit_orig.Rds")