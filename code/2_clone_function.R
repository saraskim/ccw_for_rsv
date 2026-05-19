# -----------------------------------------------------------------------------
# 
# Function for making clone dataset
#
# -----------------------------------------------------------------------------

# ----- Set up
here::i_am("code/2_clone_function.R")
pacman::p_load("dplyr","tidyr", "lubridate", "zoo", "purrr")


# ----- Function 
make_clone_dataset <- function(data, 
                               subject_id_col,
                               wt_mab_col) {
  
  data <- data %>% mutate(subject_id = .data[[subject_id_col]],
                          wt_mab = .data[[wt_mab_col]])
  
  mab_data <- data %>%
    filter(wt_mab != 0) %>%
    mutate(clone = "mab", wt = wt_mab) %>%
    arrange(subject_id, week)
  
  pref_data <- data %>%
    filter(wt_pref != 0) %>%
    mutate(clone = "pref", wt = wt_pref) %>%
    arrange(subject_id, week)
  
  control_data <- data %>%
    filter(wt_cntrl != 0) %>%
    mutate(clone = "control", wt = wt_cntrl) %>%
    arrange(subject_id, week)
  
  clone_data_stacked <- rbind(mab_data, pref_data, control_data)
  
  return(list(mab_data = mab_data, pref_data = pref_data, 
              control_data = control_data, clone_data_stacked = clone_data_stacked))
}
