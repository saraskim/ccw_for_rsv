# -----------------------------------------------------------------------------
# 
# Function for making propensity weights for:
#     Long-acting monoclonal antibodies
#     RSVpreF during pregnancy
#     Control
# 
# -----------------------------------------------------------------------------

# ----- Set up
here::i_am("code/1_weights_function.R")
pacman::p_load("dplyr","tidyr", "lubridate", "zoo", "purrr")

# ----- Function
make_weights <- function(
    data, 
    subject_id_col, 
    mab_week_col, 
    mab_tv_col, 
    grace_weeks, 
    pref_tv_col,
    mab_denom_formula, 
    pref_denom_formula, 
    ltfu_denom_formula, 
    mab_num_formula 
) {
  
  # Ensure columns read in as the observed data rather than a separate object
  data <- data %>% mutate(subject_id = .data[[subject_id_col]],
                          mab_week = .data[[mab_week_col]],
                          mab_tv = .data[[mab_tv_col]],
                          pref_tv = .data[[pref_tv_col]])
  
  # ----- Denominator Models ----- # 
  
  # -- mAb propensity model
  mab_propensity_model_data <- data %>%
    filter(pref_tv == 0, week <= mab_week) # filter to those whose moms did NOT receive RSVpreF and weeks prior and equal to when infants got mAbs
  
  # define model
  mab_propensity_model <- glm(
    mab_denom_formula,
    family = binomial(),
    data = mab_propensity_model_data
  )
  
  # predict propensity
  data <- data %>%
    mutate(mab_propensity_prediction = predict(
      mab_propensity_model, newdata = data, type = "response"
    ))
  
  
  # -- RSVpreF propensity model
  pref_propensity_model_data <- data %>% filter(week == 0) # filter to first week of life for everyone 
  
  # define model 
  pref_propensity_model <- glm(
    pref_denom_formula,
    family = binomial(),
    data = pref_propensity_model_data
  )
  
  # predict propensity
  data <- data %>%
    mutate(pref_propensity_prediction = predict(
      pref_propensity_model, newdata = data, type = "response"
    ))
  
  
  # -- LTFU (right censoring) propensity model
  ltfu_propensity_model_data <- data %>% filter(week <= ltfu_week) # this should be the same as original dataset
  
  # define model
  ltfu_propensity_model <- glm(
    ltfu_denom_formula,
    family = binomial(),
    data = ltfu_propensity_model_data
  )
  
  # predict propensity
  data <- data %>%
    mutate(ltfu_propensity_prediction = predict(
      ltfu_propensity_model, newdata = data, type = "response"
    ))
  
  # ----- Numerator Model ----- # 
  
  mab_num_propensity_model_data <- data %>% filter(mab_week <= grace_weeks) # filter to those who received mabs within grace period
  
  # define model
  mab_num_propensity_model <- glm(
    mab_num_formula,
    family = binomial(),
    data = mab_num_propensity_model_data
  )
  
  # predict propensity
  data <- data %>%
    mutate(mab_num_propensity_prediction = ifelse(
      week == grace_weeks, 
      1,
      predict(mab_num_propensity_model, newdata = data, type = "response")
    ))
  
  # ----- Make probabilities ----- # 
  
  data <- data %>%
    mutate(
      
      # denominators
      mab_propensity = case_when(
        (week < mab_week)  & (mab_week <= grace_weeks) ~ 1 - mab_propensity_prediction,
        (week == mab_week) & (mab_week <= grace_weeks) ~ mab_propensity_prediction,
        (week > mab_week)  & (mab_week <= grace_weeks) ~ 1,
        (week <= grace_weeks) & (mab_week > grace_weeks) ~ 1 - mab_propensity_prediction,
        (week > grace_weeks)  & (mab_week > grace_weeks) ~ 99999
      ),
      # mab_control_propensity = ifelse(
      #   week < mab_week, 1 - mab_propensity_prediction, 99999
      # ),
      
      # numerators
      pref_num = case_when(
        pref_tv == 0               ~ 0,
        pref_tv == 1 & mab_tv == 0 ~ 1,
        pref_tv == 1 & mab_tv == 1 ~ 0
      ),
      
      control_num = case_when(
        pref_tv == 1               ~ 0,
        pref_tv == 0 & mab_tv == 0 ~ 1,
        pref_tv == 0 & mab_tv == 1 ~ 0
      ),
     
      mab_num_propensity = case_when(
        (week < mab_week)   & (mab_week <= grace_weeks) ~ 1 - mab_num_propensity_prediction,
        (week == mab_week)  & (mab_week <= grace_weeks) ~ mab_num_propensity_prediction,
        (week > mab_week)   & (mab_week <= grace_weeks) ~ 1,
        (week < grace_weeks)  & (mab_week > grace_weeks) ~ 1 - mab_num_propensity_prediction,
        (week >= grace_weeks) & (mab_week > grace_weeks) ~ 0
      ),
      
      mab_num_ideal = ifelse(
        mab_week == 0, 1, 0
      )
    )
  
  # ----- Weights ----- #
  
  data <- data %>%
    group_by(subject_id) %>%
    arrange(week) %>%
    mutate(
      wt_mab = ((1 - pref_tv) * cumprod(mab_num_propensity)) /
        ((1 - pref_propensity_prediction) * cumprod(mab_propensity * (1 - ltfu_propensity_prediction))),
      wt_pref = pref_num /
        (pref_propensity_prediction * cumprod((1 - mab_propensity_prediction) * (1 - ltfu_propensity_prediction))),
      wt_cntrl = control_num /
        ((1 - pref_propensity_prediction) * cumprod((1 - mab_propensity_prediction) * (1 - ltfu_propensity_prediction))),
      wt_mab_ideal = (1 - pref_tv) * mab_num_ideal /
        ((1 - pref_propensity_prediction) * cumprod(mab_propensity * (1 - ltfu_propensity_prediction)))
    ) %>%
    ungroup()
  
  return(data)
}



