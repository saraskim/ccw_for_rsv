# -----------------------------------------------------------------------------
#
# Extract results
#
# -----------------------------------------------------------------------------

# --- monoclonal antibodies 

# HR
observed_estimates$msm_coef_factor.clone.mab
bootstrap_estimates %>%
  summarise(
    observed = 1 - exp(observed_estimates$msm_coef_factor.clone.mab),
    lo = 1 - exp(quantile(msm_coef_factor.clone.mab, 0.025, na.rm = TRUE)),
    hi = 1 - exp(quantile(msm_coef_factor.clone.mab, 0.975, na.rm = TRUE))
  )

# RMST
observed_estimates$rmst_mab
observed_estimates$rmst_control

bootstrap_estimates %>%
  summarise(
    observed = observed_estimates$rmst_mab,
    lo = quantile(rmst_mab, 0.025, na.rm = TRUE),
    hi = quantile(rmst_mab, 0.975, na.rm = TRUE)
  )

# plot CIs
cum_inc_plot_data <- data.frame(
  week = 1:25,
  
  # Observed estimates
  obs_control = unlist(observed_estimates[, paste0("cum_inc_control_wk", 1:25)]),
  obs_mab     = unlist(observed_estimates[, paste0("cum_inc_mab_wk",     1:25)]),
  
  # Bootstrap 95% CIs
  lo_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.025, na.rm = TRUE)),
  hi_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.975, na.rm = TRUE)),
  
  lo_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.025, na.rm = TRUE)),
  hi_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.975, na.rm = TRUE))
  
)

ggplot(cum_inc_plot_data, aes(x = week)) +
  geom_ribbon(aes(ymin = lo_control, ymax = hi_control), alpha = 0.2, fill = "gray40") +
  geom_ribbon(aes(ymin = lo_mab,     ymax = hi_mab),     alpha = 0.2, fill = "red") +
  geom_line(aes(y = obs_control, color = "Control"), lwd = 1) +
  geom_line(aes(y = obs_mab,     color = "mAb"),     lwd = 1) +
  scale_color_manual(values = c("Control" = "gray40", "mAb" = "tomato3")) +
  labs(
    x     = "Week",
    y     = "Cumulative Incidence",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "top")





# --- RSVpreF

# HR
observed_estimates$msm_coef_factor.clone.pref

bootstrap_estimates %>%
  summarise(
    observed = observed_estimates$msm_coef_factor.clone.pref,
    lo = quantile(msm_coef_factor.clone.pref, 0.025, na.rm = TRUE),
    hi = quantile(msm_coef_factor.clone.pref, 0.975, na.rm = TRUE)
  )

1-exp(-0.5456766)
1-exp(-14.92897)
1-exp(2.064847)

# RMST
observed_estimates$rmst_pref
observed_estimates$rmst_control

bootstrap_estimates %>%
  summarise(
    observed = observed_estimates$rmst_pref,
    lo = quantile(rmst_pref, 0.025, na.rm = TRUE),
    hi = quantile(rmst_pref, 0.975, na.rm = TRUE)
  )

# plot CIs
cum_inc_plot_data <- data.frame(
  week = 1:25,
  
  # Observed estimates
  obs_control = unlist(observed_estimates[, paste0("cum_inc_control_wk", 1:25)]),
  obs_pref     = unlist(observed_estimates[, paste0("cum_inc_pref_wk",     1:25)]),
  
  # Bootstrap 95% CIs
  lo_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.025, na.rm = TRUE)),
  hi_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.975, na.rm = TRUE)),
  
  lo_pref = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_pref_wk", w)]], 0.025, na.rm = TRUE)),
  hi_pref = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_pref_wk", w)]], 0.975, na.rm = TRUE))
  
)

ggplot(cum_inc_plot_data, aes(x = week)) +
  geom_ribbon(aes(ymin = lo_control, ymax = hi_control), alpha = 0.2, fill = "gray40") +
  geom_ribbon(aes(ymin = lo_pref,     ymax = hi_pref),     alpha = 0.2, fill = "cornflowerblue") +
  geom_line(aes(y = obs_control, color = "Control"), lwd = 1) +
  geom_line(aes(y = obs_pref,     color = "RSVpreF"),     lwd = 1) +
  scale_color_manual(values = c("Control" = "gray40", "RSVpreF" = "cornflowerblue")) +
  labs(
    x     = "Week",
    y     = "Cumulative Incidence",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "top")






# --- CE

# HR
obs_ratio <- exp(
  observed_estimates$msm_coef_factor.clone.mab -
    observed_estimates$msm_coef_factor.clone.pref
)

bootstrap_estimates %>%
  summarise(
    observed = obs_ratio,
    lo = quantile(exp(msm_coef_factor.clone.mab - msm_coef_factor.clone.pref), 0.025, na.rm = TRUE),
    hi = quantile(exp(msm_coef_factor.clone.mab - msm_coef_factor.clone.pref), 0.975, na.rm = TRUE)
  )

# RMST
observed_estimates$rmst_pref
observed_estimates$rmst_mab

bootstrap_estimates %>%
  summarise(
    observed = observed_estimates$rmst_pref,
    lo = quantile(rmst_pref, 0.025, na.rm = TRUE),
    hi = quantile(rmst_pref, 0.975, na.rm = TRUE)
  )

# plot CIs
cum_inc_plot_data <- data.frame(
  week = 1:25,
  
  # Observed estimates
  obs_mab = unlist(observed_estimates[, paste0("cum_inc_mab_wk", 1:25)]),
  obs_pref     = unlist(observed_estimates[, paste0("cum_inc_pref_wk",     1:25)]),
  
  # Bootstrap 95% CIs
  lo_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.025, na.rm = TRUE)),
  hi_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.975, na.rm = TRUE)),
  
  lo_pref = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_pref_wk", w)]], 0.025, na.rm = TRUE)),
  hi_pref = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_pref_wk", w)]], 0.975, na.rm = TRUE))
  
)

ggplot(cum_inc_plot_data, aes(x = week)) +
  geom_ribbon(aes(ymin = lo_mab, ymax = hi_mab), alpha = 0.2, fill = "tomato3") +
  geom_ribbon(aes(ymin = lo_pref,     ymax = hi_pref),     alpha = 0.2, fill = "cornflowerblue") +
  geom_line(aes(y = obs_mab, color = "mAb"), lwd = 1) +
  geom_line(aes(y = obs_pref,     color = "RSVpreF"),     lwd = 1) +
  scale_color_manual(values = c("mAb" = "tomato3", "RSVpreF" = "cornflowerblue")) +
  labs(
    x     = "Week",
    y     = "Cumulative Incidence",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "top")







# -- marketscan
# plot CIs
cum_inc_plot_data <- data.frame(
  week = 1:25,
  
  # Observed estimates
  obs_control = unlist(results[, paste0("cum_inc_control_wk", 1:25)]),
  obs_pref     = unlist(results[, paste0("cum_inc_pref_wk",     1:25)]),
  obs_mab = unlist(results[, paste0("cum_inc_mab_wk",     1:25)])
  
  # Bootstrap 95% CIs
  #lo_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.025, na.rm = TRUE)),
  #hi_control = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_control_wk", w)]], 0.975, na.rm = TRUE)),
  
  #lo_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.025, na.rm = TRUE)),
  #hi_mab = sapply(1:25, function(w) quantile(bootstrap_estimates[[paste0("cum_inc_mab_wk", w)]], 0.975, na.rm = TRUE))
  
)

ggplot(cum_inc_plot_data, aes(x = week)) +
  #geom_ribbon(aes(ymin = lo_control, ymax = hi_control), alpha = 0.2, fill = "gray40") +
  #geom_ribbon(aes(ymin = lo_mab,     ymax = hi_mab),     alpha = 0.2, fill = "red") +
  #geom_line(aes(y = obs_control, color = "Control"), lwd = 1) +
  geom_line(aes(y = obs_mab, color = "mAb"), lwd = 1) +
  geom_line(aes(y = obs_pref,     color = "RSVpreF"),     lwd = 1) +
  scale_color_manual(values = c("mAb" = "tomato3", "RSVpreF" = "cornflowerblue")) +
  labs(
    x     = "Week",
    y     = "Cumulative Incidence",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "top")


(1-exp(-1.620586))*100

