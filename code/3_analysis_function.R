# -----------------------------------------------------------------------------
# 
# Function for getting results:
#      Summary hazard ratios
#      KM style cumulative incidence curves over time
#      Restricted mean survival time 
#
# -----------------------------------------------------------------------------

# ----- Set up
here::i_am("code/3_analysis_function.R")
pacman::p_load("dplyr","tidyr", "lubridate", "zoo", "purrr")


# ----- function 

get_results <- function(clone_data, mab_data, control_data, pref_data, 
                        msm_formula,
                        km_formula,
                        n_weeks = 6*4.3,
                        n_wks_rmst = 6*4.3) {
  
  # -- Fit Models 
  
  msm_fit <- glm(
    msm_formula,
    weights = wt,
    data = clone_data,
    family = binomial()
  )
  
  km_control <- glm(
    km_formula,
    weights = wt_cntrl,
    data = control_data
  )
  
  km_mab <- glm(
    km_formula,
    weights = wt_mab,
    data = mab_data
  )
  
  km_pref <- glm(
    km_formula,
    weights = wt_pref,
    data = pref_data
  )
  
  # -- Predictions 
  
  pred_df <- data.frame(week = 1:n_weeks)
  
  # KM-style hazards
  haz_control <- predict(km_control, newdata = pred_df, type = "response")
  haz_mab     <- predict(km_mab,     newdata = pred_df, type = "response")
  haz_pref    <- predict(km_pref,    newdata = pred_df, type = "response")
  
  # Cumulative incidence
  cum_inc_control <- 1 - cumprod(1 - haz_control)
  cum_inc_mab     <- 1 - cumprod(1 - haz_mab)
  cum_inc_pref    <- 1 - cumprod(1 - haz_pref)
  
  # Survival
  S_control <- c(1, cumprod(1 - haz_control))
  S_mab     <- c(1, cumprod(1 - haz_mab))
  S_pref    <- c(1, cumprod(1 - haz_pref))
  
  # RMST
  rmst_control <- sum(S_control[1:n_wks_rmst])
  rmst_mab     <- sum(S_mab[1:n_wks_rmst])
  rmst_pref    <- sum(S_pref[1:n_wks_rmst])
  
  # -- Wide Results Dataframe
  
  result <- data.frame(
    # RMST
    rmst_control = rmst_control,
    rmst_mab     = rmst_mab,
    rmst_pref    = rmst_pref,
    # RMST contrasts
    rmst_mab_vs_control  = rmst_mab  - rmst_control,
    rmst_pref_vs_control = rmst_pref - rmst_control,
    rmst_mab_vs_pref     = rmst_mab  - rmst_pref,
    # Cumulative incidence at each week (KM-style)
    setNames(as.list(cum_inc_control), paste0("cum_inc_control_wk", 1:n_weeks)),
    setNames(as.list(cum_inc_mab),     paste0("cum_inc_mab_wk",     1:n_weeks)),
    setNames(as.list(cum_inc_pref),    paste0("cum_inc_pref_wk",    1:n_weeks)),
    # MSM coefficients
    setNames(as.list(coef(msm_fit)), paste0("msm_coef_", names(coef(msm_fit))))
  )
  
  return(result)
  
}

