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
  - mAb receipt under an ideal scenario
  - RSVpreF maternal vaccination
  - lost to follow-up

- To use the function, define:
  - `data`: Dataset name
  - `subject_id_col`: Variable for unique participant identifier
  - `mab_week_col`: Week mAb intervention first occurred
  - `mab_tv_col`: Time-varying mAb receipt indicator (`0 = no`, `1 = yes`)
  - `grace_weeks`: Number of weeks where not receiving intervention yet is consistent with protocol
  - `pref_tv_col`: Indicator whether subject's mother received maternal RSVpreF
  - `mab_denom_formula`: Denominator logistic regression model for time-varying mAb receipt
  - `pref_denom_formula`: Denominator logistic regression model for maternal RSVpreF receipt
  - `ltfu_denom_formula`: Denominator logistic regression model for time-varying loss to follow-up
  - `mab_num_formula`: Numerator logistic regression model for counterfactual mAb strategy

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