tbl_duration_unvacc %>% 
  filter(model == "simplified") %>% 
  ggplot(aes(t, 1-F)) +
  stat_lineribbon(aes(fill = model), .width = 0.95, size = 0.5, alpha = 0.4) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Days since infection",
    y = "Proportion still positive"
  )
ggsave("2023-03_BayesComp_ATACCC.png", width = 7, height = 7)

tbl_duration_unvacc %>% 
  filter(model == "simplified") %>% 
  group_by(.draw) %>% 
  arrange(t, .by_group = TRUE) %>% 
  mutate(S = 1 - lag(F)) %>% 
  group_by(t) %>% 
  median_qi(S, F) %>% 
  rename(time = t) %>% 
  readr::write_csv("ATACCC_posterior_summary.csv")

tbl_data = readRDS(file=here::here("Data/trajectories.RDS")) %>%
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

yrep = function(fit, fixed_fp_ct) {
  if (fixed_fp_ct) {
    fp_ct_mean = -100
    fp_ct_sd = 0.0001
    raw_draws = spread_draws(
      fit,
      pred_vl[m],
      v_s, l_fp,
      ndraws = 1000
    )
  } else {
    raw_draws = spread_draws(
      fit,
      pred_vl[m],
      v_s, l_fp,
      fp_ct_mean, fp_ct_sd, 
      ndraws = 1000
    )
  }
  raw_draws %>% 
    mutate(
      i = m %/% 42 + 1,
      day = m %% 42 - 15,
    ) %>% 
    semi_join(tbl_data, by = c("i", "day")) %>% 
    mutate(
      measure_vl = if_else(
        log(runif(n())) < -l_fp,
        rnorm(n(), fp_ct_mean, fp_ct_sd),
        rnorm(n(), pred_vl, v_s),
      ),
      measure_vl = pmax(measure_vl, 0),
    )
}

yrep(fits$simplified, TRUE) %>% 
  ggplot(aes(measure_vl)) +
  geom_density(aes(colour = "yrep", group = .draw), alpha = 0.001) +
  geom_density(aes(y, colour = "y"), data = tbl_data, size = 1) +
  labs(
    x = "log v(t)",
    colour = ""
  )

fits$simplified %>% 
  spread_draws(pred_vl[m], v_s, l_fp) %>% 
  ungroup() %>% 
  mutate(
    i = m %/% 42 + 1,
    day = m %% 42 - 15,
    yrep = rnorm(n(), pred_vl, v_s),
    # test_result = if_else(
    #   log(runif(n())) > -l_fp,
    #   rnorm(n(), pred_vl, v_s) > 0,
    #   FALSE
    # ),
  ) %>%
  filter(i <= 2) %>% 
  group_by(i, day) %>% 
  median_qi(yrep) %>% 
  left_join(tbl_data, by = c("i", "day")) %>% 
  ggplot() +
  # geom_histogram(n_pos, aes(colour = model)) +
  # geom_vline(aes(xintercept = ))
  geom_lineribbon(aes(day, yrep, ymin = .lower, ymax = .upper), alpha = 0.5, size = 0.7) +
  geom_point(aes(day, y, colour = y > 0))  +
  coord_cartesian(c(-1, 20), c(0, 20)) +
  facet_wrap(~i, ncol = 5) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Time",
    y = "Log viral load"
  )
