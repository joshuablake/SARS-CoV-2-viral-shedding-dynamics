library(bayesplot)
library(dplyr)
library(loo)
library(rstan)
library(tidybayes)
library(tidyr)
data = trajectories %>%
  mutate(group = case_when(
    (vaccinated == FALSE) & (WGS == "Pre-Alpha") ~ 1,
    (vaccinated == FALSE) & (WGS == "Alpha") ~ 1,
    (vaccinated == FALSE) & (WGS == "Delta") ~ 1,
    (vaccinated == TRUE) & (WGS == "Delta") ~ 2
  ))  %>%
  filter(!is.na(group)) %>% 
  filter(!is.na(copy)) %>% 
  mutate(
    i = as.integer(factor(id_sub, levels=unique(id_sub))),
    y = if_else(
      copy == 1,
      0,
      log(copy) - 3.43
    )
  )
yrep = spread_draws(
  fit,
  pred_vl[m],
  v_s, fp_ct_mean, fp_ct_sd, fp,
  ndraws = 1000
) %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
  ) %>% 
  semi_join(data) %>% 
  mutate(
    measure_vl = if_else(
      runif(n()) < fp,
      rnorm(n(), fp_ct_mean, fp_ct_sd),
      rnorm(n(), pred_vl, v_s),
    ),
    measure_vl = pmax(measure_vl, 0),
  ) %>% 
  pivot_wider(id_cols = .draw, values_from = measure_vl, names_from = m) %>% 
  select(-.draw)
yrep2 = spread_draws(
  fit2,
  pred_vl[m],
  v_s, l_fp,
  ndraws = 1000
) %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
  ) %>% 
  semi_join(data) %>% 
  mutate(
    measure_vl = if_else(
      log(runif(n())) < -l_fp,
      0,
      rnorm(n(), pred_vl, v_s),
    ),
    measure_vl = pmax(measure_vl, 0),
  ) %>% 
  pivot_wider(id_cols = .draw, values_from = measure_vl, names_from = m) %>% 
  select(-.draw)
ppc_dens_overlay(data$y, as.matrix(yrep))
ppc_dens_overlay(data$y, as.matrix(yrep2))
ppc_dens_overlay(data$y, as.matrix(yrep)) + xlim(0, 10)
# There's a spike in copy values around 4 which the model is struggling to explain
ppc_dens_overlay_grouped(data$y, as.matrix(yrep), data$y <= 5)
ppc_dens_overlay_grouped(data$y, as.matrix(yrep2), data$y <= 5)
ppc_stat(data$y, as.matrix(yrep), function(x) mean(x <= 0.0), binwidth = 0.001)
ppc_stat(data$y, as.matrix(yrep2), function(x) mean(x <= 0.0), binwidth = 0.001)
ppc_stat_grouped(data$y, as.matrix(yrep), data$id_sub,
                 function(x) mean(x <= 0.0), binwidth = 0.001)
ppc_stat(data$y, as.matrix(yrep2), data$id_sub,
         function(x) mean(x <= 0.0), binwidth = 0.001)
ppc_stat(data$y, as.matrix(yrep), function(x) mean(x <= 8.0), binwidth = 0.001)

draws = spread_draws(
  fit,
  pred_vl[m],
  ndraws = 1
)

draws %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
  ) %>% 
  semi_join(data)

data %>% 
  anti_join(
    draws %>% 
      mutate(
        i = m %/% 42 + 1,
        day = m %% 42 - 15,
      ) %>% 
      semi_join(data)
  )

# Outliers: test positive but shouldn't be!
spread_draws(
  fit,
  pred_vl[m],
  ndraws = 500
) %>% 
  summarise(p_pos = mean(pred_vl > 0)) %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
  ) %>% 
  right_join(data) %>% 
  filter(y > 0 & p_pos < 0.05)

# These really look like fitting issues where false negatives are being given too much weight
data %>% 
  filter(id_sub %in% c("ATA0029","ATA0039",
                       "ATA0055", "ATA0342")) %>% 
  ggplot(aes(day, copy)) +
  geom_point(aes(colour = copy > 1)) +
  scale_y_log10() +
  geom_line() +
  facet_wrap(~id_sub)

loo1 = loo(fit, save_psis = TRUE)
loo2 = loo(fit2, save_psis = TRUE)
loo_compare(loo1, loo2)

results2 = spread_draws(
  fit2,
  pred_vl[m],
  v_s, l_fp,
  ndraws = 1000
) %>% 
  ungroup() %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
  ) %>% 
  right_join(data) %>% 
  filter(y > 0) %>% 
  mutate(
    measure_vl_true = rnorm(n(), pred_vl, v_s),
    measure_vl = if_else(
      log(runif(n())) < -l_fp,
      0,
      measure_vl_true,
    ),
      yrep = pmax(measure_vl, 0),
    # Mean y conditional on a positive test
    # Formula: https://en.wikipedia.org/wiki/Folded_normal_distribution
    pred_y_pos = v_s * sqrt(2/pi) * exp(-pred_vl^2/(2*v_s^2)) -
      pred_vl * erf(-pred_vl/v_s),
    residual = y - pred_y_pos,
  )

results2 %>% 
  ggplot(aes(y, pred_y_pos)) +
  geom_point()

results2 %>% 
  ggplot(aes(y, residual)) +
  geom_point()

results2 %>% 
  group_by(i, day) %>% 
  summarise(big = mean(abs(residual) >= 5)) %>% 
  arrange(desc(big))

results2 %>% 
  filter(.draw %in% sample(unique(.draw), 12)) %>% 
  ggplot(aes(sample = residual)) +
  stat_qq_line() +
  stat_qq() +
  facet_wrap(~.draw)

results2 %>% 
  filter(.draw %in% sample(unique(.draw), 12)) %>% 
  ggplot(aes(residual)) +
  geom_histogram(aes(y = ..density..), bins = 40) +
  geom_density() +
  facet_wrap(~.draw) +
  geom_vline(xintercept = 0)

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
results2 %>% 
  group_by(.draw) %>% 
  summarise(mean = mean(residual), skew = moments::skewness(residual),
            median = median(residual)) %>% 
  pivot_longer(-.draw) %>% 
  ggplot(aes(value)) +
  geom_density() +
  geom_histogram(bins = 40) +
  facet_wrap(~name)

results2 %>% 
  filter(i %in% c(9, 49, 44, 7, 12)) %>% 
  group_by(day, i, y) %>% 
  median_qi(measure_vl_true) %>% 
  ggplot(aes(day)) +
  geom_point(aes(y = y, colour = y > 0)) +
  geom_line(aes(y = measure_vl_true)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.4) +
  facet_wrap(~i)
            