library(dplyr)
library(rstan)

prior_only = FALSE

fp_dist = FALSE

if (!prior_only) {
  df_merged <- readRDS(file="Data/trajectories.RDS") 
  df_merged <- df_merged %>%
    mutate(group = case_when(
      (vaccinated == FALSE) & (WGS == "Pre-Alpha") ~ 1,
      (vaccinated == FALSE) & (WGS == "Alpha") ~ 1,
      (vaccinated == FALSE) & (WGS == "Delta") ~ 1,
      (vaccinated == TRUE) & (WGS == "Delta") ~ 2
    ))  %>%
    filter(!is.na(group)) # Remove participants which don't fit into any group

  df_merged$obs_id <- as.integer(factor(df_merged$id_sub, levels=unique(df_merged$id_sub)))
}

# Copy number
options(mc.cores = 8)
rstan_options(auto_write = TRUE)
if (fp_dist) {
  model_file = "model.stan"
} else {
  model_file = "model2.stan"
}
model <- stan_model(model_file)

if (prior_only) {
  N = 0
  M = 0
  t = numeric()
  x = numeric()
  obs_id = integer()
  NG = 1
  group_index = integer()
} else {
  x <- df_merged$copy
  x[x==1] <- exp(3.43)
  x <- log(x) - 3.43
  t <- df_merged$day[!is.na(x)]
  obs_id <- df_merged$obs_id[!is.na(x)]
  x <- x[!is.na(x)]
  N <- length(unique(obs_id))
  group_index <- df_merged[!is.na(x),] %>%
    group_by(obs_id) %>%
    summarise(group_index= unique(group))
  group_index <- group_index$group_index
  # group_index <- rep(1, N)
  NG <- length(unique(group_index))
  M <- length(x)
}

# Pass data to model
data.stan<-list(M=M,
               N=N,
               NG=NG,
               group_index=group_index, 
               t=t,
               x=x,
               obs_id=obs_id,
               pr_fp=3,         # mean of prior for error proportion
               pr_fp_sd=0.5,      # SD of prior of error proportion
               hp_v = c(15.0, 1.25, 0.5),
               hp_v_sd = c(15.0, 0.75, 1.4),
               hp_v_sd_sd = c(10.0, 1.0, 1.0),
               pr_v_s=3.0,
               pr_v_s_sd=3.0,
               t_sd=5,
               t_mu=5,
               lod=3.4,
               eta=1
)

if (fp_dist) {
  # prior for mean of error CT distribution
  data.stan = c(
    data.stan,
    list(
      pr_fp_v=3,       
      pr_fp_v_s=3
    )
  )
}

initf <- function(chain_id) {
 list(
   v_s = 3 + (runif(1) - 0.5),
   vgrow = array(c(15, 1, 0.5), dim = c(3, NG)),
   v_sd = c(2, 0.5, 0.5), fp_ct_sd=1.51+(runif(1)-0.5)
 )
}

fit = sampling(model, chains=8, cores=8, data=data.stan,
           iter=8000, init=initf, init_r=0.1, warmup=3000, seed=497101, thin=2,
           control=list(max_treedepth=12, adapt_delta = 0.994, stepsize=0.01))

if (fp_dist) {
  saveRDS(fit, file = "fit.Rds")
} else {
  saveRDS(fit, file = "fit.Rds")
}