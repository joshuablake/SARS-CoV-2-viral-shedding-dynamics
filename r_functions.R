# Functions
copy_to_pfu_index <- function(index){
  # Convert PCR participant index to PFU participant index
  id <- unique(df_merged[df_merged$obs_id==index, "id_sub"])
  return(unique(df_pfu_id[df_pfu_id$id_sub==id, "obs_id"]))
}

y_label <- function(ii, n_cols){
  # Remove y-labels not on leftmost column
  if ((ii-1)%%(n_cols)==0){
    r <- element_text()
  } else {
    r <- element_blank()
  }
  return(r)
}

x_label <- function(ii, n_cols){
  # Remove x-labels not on bottom row
  if (ii %in% seq(N-n_cols+1, N)){
    r <- element_text()
  } else {
    r <- element_blank()
  }
  return(r)
}

posterior_pfu_i <- function(ii) {
  # Returns PFU posterior of participant ii
  ii <- copy_to_pfu_index(ii)
  return(
    data.frame(
      "v_a" = posterior_pfu$v_a[, ii],
      "v_b" = posterior_pfu$v_b[, ii],
      "p_t_max" = posterior_pfu$p_t_max[, ii],
      "p_ln_v_max" = exp(posterior_pfu$p_ln_v_max[, ii]),
      "v_s" = posterior_pfu$v_s,
      "fp" = posterior_pfu$fp
    )
  )
}

posterior_copy_i <- function(ii) {
  # Returns PCR posterior of participant ii
  return(
    data.frame(
      "v_a" = posterior_copy$v_a[, ii],
      "v_b" = posterior_copy$v_b[, ii],
      "p_t_max" = posterior_copy$p_t_max[, ii],
      "p_ln_v_max" = exp(posterior_copy$p_ln_v_max[, ii]),
      "v_s" = posterior_copy$v_s,
      "fp" = posterior_copy$fp
    )
  )
}

extract_exclude_pfu2 <- function(part) {
  # Remove trajectories to be excluded
  if (part %in% exclude_pfu) {
    mu <-
      data.frame(
        "t" = NA,
        "vl" = NA,
        "CI_l" = NA,
        "CI_u" = NA
      )
  } else {
    mu <- extract_pred_pfu(copy_to_pfu_index(part))
  }
  return(mu)
}

group_quantiles <- function(group, lod, post, measure){
  t_len_comb <- group_mean(group, lod, post, measure)
  return(data.frame("median"=quantile(t_len_comb[[1]], probs = 0.5),
                    "l_CI"= quantile(t_len_comb[[1]], probs = 0.025),
                    "u_CI"=quantile(t_len_comb[[1]], probs = 0.975)))
}

format_quantiles_2 <- function(v) {
  return(paste(sprintf("%.2f", v$median),
               " (",
               sprintf("%.2f", v$l_CI),
               ", ",
               sprintf("%.2f", v$u_CI),
               ")",
               sep = ""))
}

print_results <- function(group, group_name, lod, post) {
  return(
    data.frame(
      "group" = group_name,
      "N" = length(group),
      "t_to_peak" = format_quantiles_2(
        group_quantiles(group, lod, post, "t_to_peak")
      ),
      "t_from_peak" = format_quantiles_2(
        group_quantiles(group, lod, post, "t_from_peak")
      ),
      "v_a" = format_quantiles_2(
        group_quantiles_rates(group, post, "v_a")
      ),
      "v_b" = format_quantiles_2(
        group_quantiles_rates(group, post, "v_b")
      ),
      "peak" = format_quantiles_2(
        group_quantiles_peak(group, post, lod)
      )
    )
  )
}

group_mean <- function(group, lod, post, measure) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    time_meas <- extract_t(post(i), lod)[measure][[1]]
    df_t <- rbind(df_t, data.frame("x" = time_meas))
  }
  return(df_t)
}

group_mean_peak <- function(group, post) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    peak_val <- log10(post(i)$p_ln_v_max)
    df_t <- rbind(df_t, data.frame("x" = peak_val))
  }
  return(df_t)
}

group_mean_rates <- function(group, post, measure) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    peak_val <- post(i)[[measure]]
    df_t <- rbind(df_t, data.frame("x" = peak_val))
  }
  return(df_t)
}

group_quantiles_peak <- function(group, post, lod){
  peak_comb <- group_mean_peak(group, post)[[1]]+lod/log(10)
  return(data.frame("median"=quantile(peak_comb, probs = 0.5),
                    "l_CI"= quantile(peak_comb, probs = 0.025),
                    "u_CI"=quantile(peak_comb, probs = 0.975)))
}

group_quantiles_rates <- function(group, post, measure){
  peak_comb <- group_mean_rates(group, post, measure)[[1]]
  return(data.frame("median"=quantile(peak_comb, probs = 0.5),
                    "l_CI"= quantile(peak_comb, probs = 0.025),
                    "u_CI"=quantile(peak_comb, probs = 0.975)))
}

calc_auc <- function(part, pred) {
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:t_len)] <-
    pred[, (1:t_len) + (part - 1) * t_len]
  return(rowSums(pred_vl_df)/log(10))
}

group_auc_pfu <- function(group) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    auc <- calc_auc(copy_to_pfu_index(i), pred_pfu_samples)
    df_t <- rbind(df_t, data.frame("x" = auc))
  }
  return(df_t)
}

group_auc <- function(group) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    auc <- calc_auc(i, pred_vl_samples)
    df_t <- rbind(df_t, data.frame("x" = auc))
  }
  return(df_t)
}

auc_quantiles <- function(v){
  return(data.frame(
    "median" = quantile(v, prob=0.5),
    "l_CI" = quantile(v, prob=0.025),
    "u_CI" = quantile(v, prob=0.975)
  ))
}

print_auc_pfu <- function(group, group_name) {
  v <- group_auc_pfu(group)[[1]]
  return(data.frame(
    "group" = group_name,
    "N" = length(group),
    "AUC" = format_quantiles_2(auc_quantiles(v))
  ))
}

print_auc_copy <- function(group, group_name) {
  v <- group_auc(group)[[1]]
  return(data.frame(
    "group" = group_name,
    "N" = length(group),
    "AUC" = format_quantiles_2(auc_quantiles(v))
  ))
}

extract_t <- function(posterior, lod) {
  p_t_max <- posterior$p_t_max
  v_a <- posterior$v_a
  v_b <- posterior$v_b
  p_ln_v_max <- log(posterior$p_ln_v_max)
  t_start <-
    p_t_max - (log(v_a + v_b) - log(v_b) + p_ln_v_max) / v_a
  t_end <-
    p_t_max + (log(v_a + v_b) - log(v_a) + p_ln_v_max) / v_b
  t_to_peak <- p_t_max - t_start
  t_from_peak <- t_end - p_t_max
  t_length <- t_end - t_start
  return(
    data.frame(
      "t_start" = t_start,
      "t_end" = t_end,
      "t_to_peak" = t_to_peak,
      "t_from_peak" = t_from_peak,
      "t_length" = t_length
    )
  )
}

calc_pp <- function(p1, p2, x_lim) {
  p1_den <-
    density(p1,
            from = -x_lim,
            to = x_lim)
  p2_den <-
    density(p2,
            from = -x_lim,
            to = x_lim)
  mdx <- mean(diff(p1_den$x))
  rx <-
    seq(-mdx * (length(p1_den$x) - 1), mdx * (length(p1_den$x) -
                                                1), by = mdx)
  r <-
    convolve(p2_den$y, # longer
             p1_den$y, # shorter
             conj = TRUE,
             type = "open")
  pp <-
    data.frame("x" = rx,
               "y" = r)
  pp$y <- pp$y / integrate.xy(pp$x, pp$y)
  return(pp)
}

calc_conv <- function(den1, den2) {
  mdx <- mean(diff(den1$x))
  rx <-
    seq(-mdx * (length(den1$x) - 1), mdx * (length(den1$x) -
                                              1), by = mdx)
  r <-
    convolve(den2$y, # longer
             den1$y, # shorter
             conj = TRUE,
             type = "open")
  pp <-
    data.frame("x" = rx,
               "y" = r)
  pp$y <- pp$y / integrate.xy(pp$x, pp$y)
  return(pp)
}

pp_int <- function(pp){
  return(integrate.xy(pp$x[pp$x>0], pp$y[pp$x>0]))
}

bayes_fac <- function(p){
  p_int <- pp_int(p)
  return(p_int/(1-p_int))
}


##

extract_pred <- function(part) {
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:t_len)] <-
    pred_vl_samples[, (1:t_len) + (part - 1) * t_len] + 3.4
  
  pred_vl_0.5 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  pred_vl_0.05 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  pred_vl_0.95 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = tr,
      "vl" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

extract_data <- function(part){
  df_part <- df_merged[df_merged$obs_id==part, ]
  df_part$copy[df_part$copy==1] <- exp(3.4)
  t <- df_part$day
  vl <- log(df_part$copy)
  return(data.frame("t"=t, "vl"=vl))
}

y_label <- function(ii, n_cols){
  if ((ii-1)%%(n_cols)==0){
    r <- element_text()
  } else {
    r <- element_blank()
  }
  return(r)
}

x_label <- function(ii, n_cols){
  if (ii %in% seq(N-n_cols+1, N)){
    r <- element_text()
  } else {
    r <- element_blank()
  }
  return(r)
}

extract_pred <- function(part) {
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:t_len)] <-
    pred_vl_samples[, (1:t_len) + (part - 1) * t_len] + 3.4
  
  pred_vl_0.5 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  pred_vl_0.05 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  pred_vl_0.95 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = tr,
      "vl" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

# extract_data <- function(part){
#   df_part <- df_merged[df_merged$obs_id==part, ]
#   df_part$copy[df_part$copy==1] <- exp(3.4)
#   t <- df_part$day
#   vl <- log(df_part$copy)
#   return(data.frame("t"=t, "vl"=vl))
# }

extract_pred_pfu <- function(part) {
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:t_len)] <-
    pred_pfu_samples[, (1:t_len) + (part - 1) * t_len] + 2.3
  
  pred_vl_0.5 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  pred_vl_0.05 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  pred_vl_0.95 <- pred_vl_df %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = tr,
      "vl" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

extract_data_pfu <- function(part){
  df_part <- df_pfu[df_pfu$obs_id==part, ]
  df_part$pfu[df_part$pfu==10] <- 10.01
  df_part$pfu[df_part$pfu==1] <- 10.0
  t <- df_part$day
  pfu <- log(df_part$pfu)
  return(data.frame("t"=t, "vl"=pfu))
}

extract_exclude_data_pfu <- function(part) {
  if (part %in% exclude_pfu) {
    r <- data.frame("t" = NA, "vl" = NA)
  } else{
    df_part <- df_pfu[df_pfu$obs_id == part,]
    df_part$pfu[df_part$pfu==10] <- 10.01
    df_part$pfu[df_part$pfu==1] <- 10.0
    t <- df_part$day
    pfu <- log(df_part$pfu)
    r <- data.frame("t" = t, "vl" = pfu)
  }
  return(r)
}

plot_auc_copy <- function(gplot, group, col) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    auc <- calc_auc(i, pred_vl_samples)
    df_t <- rbind(df_t, data.frame("x" = auc))
    gplot <-
      gplot + geom_density(
        data = data.frame("x" = auc),
        aes(x = x),
        fill = col,
        colour = NA,
        alpha = 0.25
      )
  }
  gplot <-
    gplot + geom_density(data = data.frame("x" = df_t), aes(x = x), colour = "black")
  return(gplot)
}

plot_auc_pfu <- function(gplot, group, col) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    auc <- calc_auc(copy_to_pfu_index(i), pred_pfu_samples)
    df_t <- rbind(df_t, data.frame("x" = auc))
    gplot <-
      gplot + geom_density(
        data = data.frame("x" = auc),
        aes(x = x),
        fill = col,
        colour = NA,
        alpha = 0.25
      )
  }
  gplot <-
    gplot + geom_density(data = data.frame("x" = df_t), aes(x = x), colour = "black")
  return(gplot)
}

calc_auc <- function(part) {
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:t_len)] <-
    pred_vl_samples[, (1:t_len) + (part - 1) * t_len]
  return(rowSums(pred_vl_df))
}

calc_corr <- function(i, j, post, positive) {
  c_vals <- post$Lc[, i, j]
  den <- density(c_vals, from=-2, to=2)
  den$y <- den$y / integrate.xy(den$x, den$y)
  if (positive){
    pp <- integrate.xy(den$x[den$x>0], den$y[den$x>0])
  } else {
    pp <- integrate.xy(den$x[den$x<0], den$y[den$x<0])
  }
  return(
    data.frame(
      "entry" = paste0(i, j),
      "med" = quantile(c_vals, prob = 0.5),
      "l_CI" = quantile(c_vals, prob = 0.025),
      "u_CI" = quantile(c_vals, prob = 0.975),
      "BF" = pp/(1-pp)
    )
  )
}

corr_bf <- function(i,j, post, positive){
  c_vals <- post$Lc[, i, j]
  den <- density(c_vals, from=-2, to=2)
  den$y <- den$y / integrate.xy(den$x, den$y)
  if (positive){
    pp <- integrate.xy(den$x[den$x>0], den$y[den$x>0])
  } else {
    pp <- integrate.xy(den$x[den$x<0], den$y[den$x<0])
  }
  return(pp)
}

compare_kinetics <- function(group1, group2, lod, post, x_lim, order) {
  measures <- c("t_to_peak", "t_from_peak", "t_length")
  c <- 0
  for (meas_i in measures) {
    c <- c + 1
    if (order[c] == 2){
      group1_c <- group2
      group2_c <- group1
    } else {
      group1_c <- group1
      group2_c <- group2
    }
    p1 <- group_mean(group1_c, lod, post, meas_i)
    p2 <-
      group_mean(group2_c,
                 lod,
                 post,
                 meas_i)
    pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
    if (c == 1) {
      res <- data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12))
    } else {
      res <- rbind(res, data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12)))
    }
  }
  return(res)
}

compare_rates <- function(group1, group2, post, x_lim) {
  measures <- c("v_a", "v_b") 
  c <- 0
  for (meas_i in measures) {
    c <- c+1
    p1 <- group_mean_rates(group1, post, meas_i)
    p2 <-
      group_mean_rates(group2,
                       post,
                       meas_i)
    
    pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
    if (c == 1) {
      res <- data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12))
    } else {
      res <- rbind(res, data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12)))
    }
  }
  return(res)
}

compare_copy_pfu_rates <- function(group1, group2, x_lim) {
  measures <- c("v_a", "v_b")
  c <- 0
  for (meas_i in measures) {
    c <- c + 1
    p1 <- group_mean_rates(group1, posterior_copy_i, meas_i)
    p2 <-
      group_mean_rates(group2,
                       posterior_pfu_i,
                       meas_i)
    pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
    if (c == 1) {
      res <- data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12))
    } else {
      res <- rbind(res, data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12)))
    }
  }
  return(res)
}

compare_peaks <- function(group1, group2, post, x_lim) {
  p1 <- group_mean_peak(group1, post)
  p2 <-
    group_mean_peak(group2,
                    post)
  pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
  
  res <-
    data.frame("measure" = "peaks", "bayes_factor" = bayes_fac(pp12))
  
  return(res)
}

compare_copy_pfu_kinetics <- function(group1, group2, x_lim) {
  measures <- c("t_to_peak", "t_from_peak", "t_length")
  c <- 0
  for (meas_i in measures) {
    c <- c + 1
    p1 <- group_mean(group1, 2.3, posterior_pfu_i, meas_i)
    p2 <-
      group_mean(group2,
                 3.4,
                 posterior_copy_i,
                 meas_i)
    pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
    if (c == 1) {
      res <- data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12))
    } else {
      res <- rbind(res, data.frame("measure"=meas_i, "bayes_factor" = bayes_fac(pp12)))
    }
  }
  return(res)
}

compare_auc_pfu <- function(group1, group2, x_lim) {
  p1 <- group_auc_pfu(group1)
  p2 <-
    group_auc_pfu(group2)
  pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
  res <-
    data.frame("measure" = "AUC PFU", "bayes_factor" = bayes_fac(pp12))
  return(res)
}

compare_auc_copy <- function(group1, group2, x_lim) {
  # group2 - higher auc
  p1 <- group_auc(group1)
  p2 <-
    group_auc(group2)
  pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
  res <-
    data.frame("measure" = "AUC Copy", "bayes_factor" = bayes_fac(pp12))
  return(res)
}

compare_auc_copy_pfu <- function(group1, group2, x_lim) {
  # group2 - higher auc
  p1 <- group_auc_pfu(group1)
  p2 <-
    group_auc(group2)
  pp12 <- calc_pp(p1$x, p2$x, x_lim = x_lim)
  res <-
    data.frame("measure" = "AUC Copy", "bayes_factor" = bayes_fac(pp12))
  return(res)
}

group_plot <- function(gplot, group, lod, post, measure, col) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    # time_meas <- extract_t(post(i), lod)[measure][[1]]
    time_meas <- post(i)[[measure]]
    df_t <- rbind(df_t, data.frame("x" = time_meas))
    gplot <-
      gplot + geom_density(
        data = data.frame("x" = time_meas),
        aes(x = x),
        fill = col,
        colour = NA,
        alpha = 0.25
      )
  }
  gplot <-
    gplot + geom_density(data = data.frame("x" = df_t), aes(x = x), colour = "black")
  return(gplot)
}

group_plot_duration <- function(gplot, group, lod, post, measure, col) {
  df_t <- data.frame("x" = double())
  for (i in group) {
    time_meas <- extract_t(post(i), lod)[measure][[1]]
    df_t <- rbind(df_t, data.frame("x" = time_meas))
    gplot <-
      gplot + geom_density(
        data = data.frame("x" = time_meas),
        aes(x = x),
        fill = col,
        colour = NA,
        alpha=0.25
      )
  }
  gplot <-
    gplot + geom_density(data = data.frame("x" = df_t), aes(x = x), colour = "black")
  return(gplot)
}

param_quantiles <- function(v){
  return(data.frame("median"=quantile(v, probs = 0.5),
                    "l_CI"= quantile(v, probs = 0.025),
                    "u_CI"=quantile(v, probs = 0.975)))
}
