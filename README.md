# RSV preventive intervention effectiveness using clone censor weighting

> This set of code provides functions for creating propensity weights, cloning the dataset, running analyses, and bootstrap.  

------------------------------------------------------------------------

## Data set up
- Dataset should be in long weekly format, with one row per week per infant
- Functions assume that the last row per infant is the week of the first occurrence of an outcome, lost to follow-up, or study end

- At minimum, dataset should have:
  - `subject_id`: Unique participant identifier
  - `week`: Week from time zero
  - `mab_tv`: Time-varying mAb receipt indicator (0 = no, 1 = yes)
  - `mab_week`: Week mAb intervention first occurred (same across all weeks)
  - `pref_tv`: Maternal RSVpreF receipt indicator (0 = no, 1 = yes, same across all weeks because received prior to birth)
  - `pref_week`: Week mAb intervention first occurred (if received, should be 0)
  - `outcome_tv`: Time-varying outcome indicator
  - `outcome_week`: Week outcome first occurred
  - `ltfu_tv`: Lost-to-follow-up indicator
  - `ltfu_week`: Week lost to follow-up occurred
  - `99999`: Indicates the event never occurred during follow-up
  
#### Example Dataset Format

| subject_id | week | mab_tv | mab_week | pref_tv | pref_week | outcome_tv | outcome_week | ltfu_tv | ltfu_week |
|------------|------|--------|-----------|----------|------------|------------|---------------|----------|------------|
| 1 | 0 | 0 | 2 | 0 | 99999 | 0 | 4 | 0 | 99999 |
| 1 | 1 | 0 | 2 | 0 | 99999 | 0 | 4 | 0 | 99999 |
| 1 | 2 | 1 | 2 | 0 | 99999 | 0 | 4 | 0 | 99999 |
| 1 | 3 | 1 | 2 | 0 | 99999 | 0 | 4 | 0 | 99999 |
| 1 | 4 | 1 | 2 | 0 | 99999 | 1 | 4 | 0 | 99999 |
| 2 | 0 | 0 | 99999 | 1 | 0 | 0 | 99999 | 0 | 2 |
| 2 | 1 | 0 | 99999 | 1 | 0 | 0 | 99999 | 0 | 2 |
| 2 | 2 | 0 | 99999 | 1 | 0 | 0 | 99999 | 1 | 2 |

## Libraries
- Please ensure the following libraries are installed and loaded:
  - dplyr
  - tidyr
  - lubridate
  - zoo
  - purrr


## Code

#### `code/1_weights_function.R`

- Ensure correctly formatted dataset is read into R.
- Function fits propensity models and makes weights for:
  - mAb receipt within the grace period
  - mAb receipt under an ideal scenario (infants receive mAbs in week 0)
  - RSVpreF maternal vaccination
  - lost to follow-up

- To use the function, define:
  - `data`: Dataset name
  - `subject_id_col`: Variable for unique participant identifier
  - `mab_week_col`: Week mAb intervention first occurred (same across all weeks)
  - `mab_tv_col`: Time-varying mAb receipt indicator (`0 = no`, `1 = yes`)
  - `grace_weeks`: Number of weeks where not receiving intervention yet is consistent with intervention protocol
  - `pref_tv_col`: Indicator whether subject's mother received maternal RSVpreF (same across all weeks since determined before birth)
  - `mab_denom_formula`: Denominator logistic regression model for time-varying monoclonal antibody receipt
  - `pref_denom_formula`: Denominator logistic regression model for maternal RSVpreF receipt
  - `ltfu_denom_formula`: Denominator logistic regression model for time-varying loss to follow-up
  - `mab_num_formula`: Numerator logistic regression model for counterfactual scenario where infants receive mAbs under the defined strategy

- Returns dataset with weights.

- Example use:

```r
data_weekly <- readRDS("data/exampledata.rds")

data_weekly_weights <- make_weights(
  data = data_weekly,
  subject_id_col = "infant_id",
  mab_week_col = "mab_week",
  mab_tv_col = "mab_tv",
  grace_weeks = 6 * 4.3,
  pref_tv_col = "pref_binary",

  mab_denom_formula =
    mab_tv ~ maternal_age + splines::ns(week) + season,

  pref_denom_formula =
    pref_tv ~ maternal_age + season,

  ltfu_denom_formula =
    ltfu_tv ~ splines::ns(epi_week_tv),

  mab_num_formula =
    mab_tv ~ maternal_age + splines::ns(week) + season
)
```

- Can save `data_weekly_weights` as a dataset if needed.


#### `code/2_clone_function.R`

- Ensure long dataset with weights produced in `code/1_weights_function.R` is read in.
- Function creates cloned dataset by stacking subsets with non-zero weights for:
  - mAbs
  - RSVpreF
  - lost to follow-up

- To use the function, define:
  - `data`: Dataset name
  - `subject_id_col`: Variable for unique participant identifier
  - `wt_mab_col`: Define which mAb strategy is of interest:
    - mAb within grace period
    - mAb under ideal scenario

- Returns an object with:
  - subsets for mAbs, RSVpreF, and lost to follow-up
  - a cloned dataset where these subsets are stacked into one dataset

- Example use:

```r
clone_data <- make_clone_dataset(
  data = data_weekly_weights,
  subject_id_col = "subject_id",
  wt_mab_col = "wt_mab"
)
```

- Can save `clone_data` if needed.


#### `code/3_analysis_function.R`

- Analyses obtained include:
  - coefficients from the marginal structural model
  - Kaplan-Meier style cumulative incidence curves over time
  - restricted mean survival time

- Ensure `clone_data` from `code/2_clone_function.R` is read in.

- To use the function, define:
  - `clone_data`: Fully stacked cloned data created in `make_clone_dataset()`
  - `mab_data`: Subset of non-zero mAb weights created in `make_clone_dataset()`
  - `control_data`: Subset of non-zero control weights created in `make_clone_dataset()`
  - `pref_data`: Subset of non-zero RSVpreF weights created in `make_clone_dataset()`
  - `msm_formula`: MSM formula for hazard ratio
  - `km_formula`: Formula for Kaplan-Meier style cumulative incidence curves
  - `n_weeks`: Number of weeks to follow-up for cumulative incidence curves
  - `n_wks_rmst`: Number of weeks to follow-up for restricted mean survival time

- Returns analysis results in a wide format (one row) to prepare for the bootstrap process.

- Example use:

```r
results <- get_results(
  clone_data = clone_data$clone_data_stacked,
  mab_data = clone_data$mab_data,
  control_data = clone_data$control_data,
  pref_data = clone_data$pref_data,
  msm_formula = rsv_hosp_tv ~ splines::ns(week) + factor(clone),
  km_formula = rsv_hosp_tv ~ factor(week),
  n_weeks = 6 * 4.3,
  n_wks_rmst = 6 * 4.3
)
```

- Can save `results` if needed.


#### `code/4_bootstrap.R`

- Runs bootstrap process.
- Ensure all functions from:
  - `code/1_weights_function.R`
  - `code/2_clone_function.R`
  - `code/3_analysis_function.R`

  are read in.

- Ensure correctly formatted weekly dataset is read in.

- Returns observed and bootstrapped results where each row contains results from one bootstrap iteration.

- Example use:

```r
boot_output <- run_bootstrap(
  data = data_weekly,
  n_boot = 10,
  seed = 123,
  subject_id_col = "infant_id",
  mab_week_col = "mab_week",
  mab_tv_col = "mab_tv",
  grace_weeks = 6 * 4.3,
  pref_tv_col = "pref_tv",
  wt_mab_col = "wt_mab",
  n_weeks = 6 * 4.3,
  n_wks_rmst = 6 * 4.3,

  mab_denom_formula =
    mab_tv ~ maternal_age + splines::ns(week) + season,

  pref_denom_formula =
    pref_tv ~ maternal_age + season,

  ltfu_denom_formula =
    ltfu_tv ~ splines::ns(epi_week_tv),

  mab_num_formula =
    mab_tv ~ maternal_age + splines::ns(week) + season,

  msm_formula =
    outcome_tv ~ splines::ns(week) + factor(clone),

  km_formula =
    outcome_tv ~ factor(week)
)
```

- NOTE:
  - This function will likely need to be run on a high-performance computer with additional optimization methods