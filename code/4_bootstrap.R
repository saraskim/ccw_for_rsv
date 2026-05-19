# -----------------------------------------------------------------------------
# 
# Function for bootstrap process
# 
# -----------------------------------------------------------------------------

# ----- Set up
here::i_am("code/4_bootstrap.R")
pacman::p_load("dplyr","tidyr", "lubridate", "zoo")

# ----- Function
run_bootstrap <- function(
    data,
    n_boot,
    seed,
    subject_id_col,
    mab_week_col,
    mab_tv_col,
    grace_weeks,
    pref_tv_col,
    wt_mab_col,
    n_weeks,
    n_wks_rmst,
    mab_denom_formula,
    pref_denom_formula,
    ltfu_denom_formula,
    mab_num_formula,
    msm_formula,
    km_formula
) {
  
  data <- data %>% mutate(subject_id = .data[[subject_id_col]])
  
  set.seed(seed)
  infant_ids <- unique(data$subject_id)
  
  # ----- Observed Estimates -----
  
  message("Computing observed estimates...")
  
  observed <- tryCatch({
    obs_weights <- make_weights(
      data                 = data,
      subject_id_col       = subject_id_col,
      mab_week_col         = mab_week_col,
      mab_tv_col           = mab_tv_col,
      grace_weeks          = grace_weeks,
      pref_tv_col          = pref_tv_col,
      mab_denom_formula    = mab_denom_formula,
      pref_denom_formula   = pref_denom_formula,
      ltfu_denom_formula   = ltfu_denom_formula,
      mab_num_formula      = mab_num_formula
    )
    
    
    obs_cloned <- make_clone_dataset(data = obs_weights,
                                     subject_id_col = subject_id_col,
                                     wt_mab_col = wt_mab_col)
    
    get_results(
      clone_data   = obs_cloned$clone_data_stacked,
      mab_data     = obs_cloned$mab_data,
      control_data = obs_cloned$control_data,
      pref_data    = obs_cloned$pref_data,
      msm_formula  = msm_formula,
      km_formula   = km_formula,
      n_weeks      = n_weeks,
      n_wks_rmst   = n_wks_rmst
    ) %>% mutate(iteration = 0, type = "observed")
    
  }, error = function(e) {
    message("Observed estimates failed: ", e$message)
    NULL
  })
  
  # ----- Bootstrap Iterations -----
  
  boot_results <- vector("list", n_boot)
  
  for (i in 1:n_boot) {
    
    if (i %% 100 == 0) message("Bootstrap iteration ", i, " of ", n_boot)
    
    boot_results[[i]] <- tryCatch({
      
      # Resample infants with replacement and properly duplicate rows
      boot_ids  <- sample(infant_ids, size = length(infant_ids), replace = TRUE)
      boot_data <- map_dfr(boot_ids, ~data %>% filter(subject_id == .x)) %>%
        group_by(subject_id) %>%
        mutate(subject_id = cur_group_id()) %>%
        ungroup()
      
      # Chain the three functions
      boot_weights <- make_weights(
        data                 = boot_data,
        subject_id_col       = subject_id_col,
        mab_week_col         = mab_week_col,
        mab_tv_col           = mab_tv_col,
        grace_weeks          = grace_weeks,
        pref_tv_col          = pref_tv_col,
        mab_denom_formula    = mab_denom_formula,
        pref_denom_formula   = pref_denom_formula,
        ltfu_denom_formula   = ltfu_denom_formula,
        mab_num_formula      = mab_num_formula
      )
      
      boot_cloned <- make_clone_dataset(data = boot_weights,
                                        subject_id_col = subject_id_col,
                                        wt_mab_col = wt_mab_col)
      
      get_results(
        clone_data   = boot_cloned$clone_data,
        mab_data     = boot_cloned$mab_data,
        control_data = boot_cloned$control_data,
        pref_data    = boot_cloned$pref_data,
        msm_formula  = msm_formula,
        km_formula   = km_formula,
        n_weeks      = n_weeks,
        n_wks_rmst   = n_wks_rmst
      ) %>% mutate(iteration = i, type = "bootstrap")
      
    }, error = function(e) {
      message("Iteration ", i, " failed: ", e$message)
      NULL
    })
  }
  
  # ----- Combine & Return -----
  
  all_results <- bind_rows(observed, boot_results)
  
  return(all_results)
}

# KPGA
system.time(
boot_output <- run_bootstrap (
    data = data_weekly,
    n_boot = 10,
    seed = 123,
    subject_id_col = "infant_id",
    mab_week_col = "mab_week",
    mab_tv_col = "mab_tv",
    grace_weeks = 6*4.3,
    pref_tv_col = "pref_tv",
    n_weeks = 6*4.3,
    n_wks_rmst = 6*4.3,
    mab_denom_formula = mab_tv ~ maternal_age + prenatal_encs_impute + preterm_impute +
      circulatory_malform_binary + pulmonary_mal_binary +
      splines::ns(week)*month_rel_season + season,
    
    pref_denom_formula = pref_tv ~ maternal_age + prenatal_encs_impute +
      month_rel_season + season,
    
    ltfu_denom_formula = ltfu_tv ~ maternal_age + prenatal_encs_impute + preterm_impute +
      circulatory_malform_binary + pulmonary_mal_binary +
      splines::ns(week) + splines::ns(epi_week_tv),
    
    mab_num_formula = mab_tv ~ maternal_age + prenatal_encs_impute + preterm_impute +
      circulatory_malform_binary + pulmonary_mal_binary +
      splines::ns(week)*month_rel_season + season,
    
    msm_formula = outcome_tv ~ splines::ns(week) + factor(clone),
    
    km_formula  = outcome_tv ~ factor(week)
)
)




