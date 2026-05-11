#' @importFrom magrittr %>%
#' @importFrom dplyr across filter group_by row_number select slice ungroup
#' @importFrom tidyr drop_na
#' @importFrom tidyselect where
#' @importFrom data.table := as.data.table fwrite rbindlist
#' @importFrom tibble as_tibble
#' @importFrom readr read_rds write_rds
#' @importFrom lme4 VarCorr
#' @importFrom stats binomial quantile sd var
#' @importFrom utils read.csv
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".data",
  "..dxlist",
  "..dxlist_main",
  "age",
  "age_k1",
  "age_k2",
  "BYr_k1",
  "BYr_k2",
  "cnuRel",
  "cohort_k1",
  "cohort_k2",
  "current_rows",
  "ID",
  "male_k2",
  "matID_k1",
  "matID_k2",
  "obs_num",
  "outcome",
  "outcome_vars",
  "patID_k1",
  "patID_k2",
  "range_max",
  "same_matID",
  "same_patID"
))
