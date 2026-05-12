#------------------------------------------------------------------------------
# Author Mason Garrison
# Date: 815
# Filename: functions_CorrelateOutcomesByGroup
# Purpose: this code calcuates the correlation between the outcomes for each kin group, groups by mtdna, cnu, and bins of R

#------------------------------------------------------------------------------

## First are the helper functions
options(scipen = 10, digits = 11)


#' Not-In Operator
#'
#' Returns \code{TRUE} for each element of \code{x} that is \emph{not} present
#' in \code{table}.  Equivalent to \code{!(\%in\%)}.
#'
#' @param x A vector of values to test.
#' @param table A vector of values to test against.
#'
#' @return A logical vector the same length as \code{x}.
#'
#' @keywords internal
`%notin%` <- Negate(`%in%`)

#' Convert Numeric Columns to Integer Where Possible
#'
#' Converts numeric columns of a data frame (or tibble) to integer when every
#' unique non-missing value in that column is a whole number.  Two modes of
#' operation are provided, selected by \code{memory_manage}.
#'
#' @param tbl A data frame or tibble.
#' @param memory_manage Integer flag controlling memory strategy.  When
#'   \code{< 1} a \pkg{dplyr}/\pkg{tidyr} approach is used; when \code{>= 1}
#'   a \pkg{data.table} in-place approach is used instead.  Default is
#'   \code{0L}.
#'
#' @return A tibble with integer columns wherever conversion was possible.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate across if_else
#' @importFrom tidyselect where
#' @importFrom data.table ":=" as.data.table
#' @importFrom tibble as_tibble
#' @keywords internal
convert_to_integer <- function(tbl, memory_manage = 0L) {
  if (memory_manage < 1) {
    tbl %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ dplyr::if_else(. == floor(.), as.integer(.), .)))
  } else {
    # convert to datatable
    tbl <- data.table::as.data.table(tbl)
    # convert numeric to integer where possible
    for (col in names(tbl)) {
      if (is.numeric(tbl[[col]])) {
        # if intger the same for both, then convert in place
        if (mean(floor(unique(tbl[[col]])) == unique(tbl[[col]]), na.rm = TRUE)) {
          tbl[, (col) := as.integer(tbl[[col]])]
        }
      }
    }
    tibble::as_tibble(tbl)
  }
}

#' Try an Expression, Return \code{NA} on Error
#'
#' Evaluates \code{expr} and returns \code{NA} (rather than stopping) if the
#' expression throws an error.
#'
#' @param expr An R expression to evaluate.
#'
#' @return The result of \code{expr}, or \code{NA} if an error occurs.
#'
#' @keywords internal
try_NA <- function(expr) {
  tryCatch(expr, error = function(err) NA)
}

#' Optionally Slice a Table to the First 1 000 Rows
#'
#' A thin wrapper around \code{dplyr::slice()} used to subsample large tables
#' during testing or memory-constrained runs.  When \code{slice_1000 = FALSE}
#' the input table is returned unchanged.
#'
#' @param tbl A data frame or tibble.
#' @param slice_1000 Logical.  If \code{TRUE}, only the first 1 000 rows are
#'   returned; if \code{FALSE} the entire table is returned.
#' @param memory_manage Integer memory-management flag (currently unused but
#'   kept for API consistency).  Default is \code{0L}.
#'
#' @return A data frame or tibble with at most 1 000 rows when
#'   \code{slice_1000 = TRUE}, otherwise the original \code{tbl}.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr slice
#' @keywords internal
sliceFunction <- function(tbl, slice_1000, memory_manage = 0L) {
  if (slice_1000 == TRUE) {
    tbl %>% dplyr::slice(1:1000)
  } else {
    tbl
  }
}
#' @rdname sliceFunction
#' @keywords internal
sliceFuction <- sliceFunction

# passes thru multiple grouping variables. right now I'm using groupby, but this could be replaced by expand grid.
# it does allow for multiple groups
#' Group a Table by One or More Variables
#'
#' A flexible wrapper around \code{dplyr::group_by()} that accepts a vector of
#' grouping variable names.  When no valid grouping variables are provided the
#' table is returned ungrouped.
#'
#' @param tbl A data frame or tibble.
#' @param grouping_vars A character vector of column names to group by, or
#'   \code{NA}/\code{NULL}/\code{0} to skip grouping.
#' @param verbose Logical.  If \code{TRUE}, messages describing the grouping
#'   operation are emitted.  Default is \code{FALSE}.
#' @param memory_manage Integer memory-management flag (currently unused but
#'   kept for API consistency).  Default is \code{0L}.
#'
#' @return A grouped (or ungrouped) data frame / tibble.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @keywords internal
group_byFunction <- function(tbl, grouping_vars, verbose = FALSE, memory_manage = 0L) {
  if (length(grouping_vars) > 1) {
    if (verbose) {
      message(paste0("grouping vars", grouping_vars))
    }
    grouping_call <- paste0(
      "tbl %>% group_by(",
      paste(grouping_vars, collapse = ", "), ")"
    )

    eval(parse(text = grouping_call))
  } else if (length(grouping_vars) == 0 || is.na(grouping_vars) || is.null(grouping_vars) || grouping_vars == 0) {
    if (verbose) {
      message(paste0("grouping vars null", grouping_vars))
    }
    tbl
  } else {
    if (verbose) {
      message(paste0("grouping vars", grouping_vars))
    }
    grouping_call <- paste0(
      "tbl %>% group_by(",
      paste(grouping_vars, collapse = ", "), ")"
    )

    eval(parse(text = grouping_call))
  }
}


# this function lets you implement other data cleaning
#' Apply Pre-Defined Mutation (Data Transformation) Steps
#'
#' Adds derived grouping columns to \code{tbl} based on the value of
#' \code{mutate_vars}.  Each value of \code{mutate_vars} maps to a specific
#' set of column-creation operations (e.g. gender groupings, linkage type,
#' birth-cohort groupings, case-control status).  If \code{mutate_vars} is
#' \code{NULL}, \code{NA}, or \code{0} the table is returned unchanged.
#'
#' @param tbl A data frame or tibble containing the kinship-pair data.
#' @param mutate_vars A single string (or \code{NULL}/\code{NA}/\code{0})
#'   naming the transformation to apply.  Supported values include:
#'   \describe{
#'     \item{\code{"gender_groupings"}}{Add a \code{gender_groupings} column
#'       coding same-sex and opposite-sex pairs.}
#'     \item{\code{"linkagetype"}}{Add columns for shared maternal/paternal IDs
#'       and a composite \code{linkagetype} variable.}
#'     \item{\code{"linkage_any"}}{Like \code{"linkagetype"} but codes any
#'       shared lineage as the same value.}
#'     \item{\code{"same_patID"}}{Add a \code{same_patID} indicator.}
#'     \item{\code{"same_matID"}}{Add a \code{same_matID} indicator.}
#'     \item{\code{"gender_groupings_linkagetype"} / \code{"gender_linkagetype"}}{
#'       Combine linkage type and gender groupings.}
#'     \item{\code{"casecontrol_groupings"}}{Add a case/control grouping column.}
#'     \item{\code{"cohort_groupings_19"}}{Add birth-cohort groupings (18th vs
#'       19th century).}
#'     \item{\code{"cohort_groupings_19flat"}}{Flattened 18th/19th-century cohort
#'       groupings.}
#'     \item{\code{"cohort_groupings_match"}}{Indicator for pairs in the same
#'       birth century.}
#'     \item{\code{"cohort_groupings"}}{Full century-level cohort groupings.}
#'     \item{\code{"gender_cohort_groupings_19"}}{Gender + 19th-century cohort
#'       groupings.}
#'     \item{\code{"gender_cohort_groupings"}}{Gender + full cohort groupings.}
#'     \item{\code{"cohort_linkagetype_19"}}{Linkage type + 19th-century cohort.}
#'     \item{\code{"cohort_linkagetype"}}{Linkage type + full cohort groupings.}
#'     \item{\code{"cohort_gender_linkagetype_19"}}{Linkage + cohort + gender
#'       (19th-century).}
#'     \item{\code{"cohort_gender_linkagetype"}}{Linkage + cohort + gender (full).}
#'     \item{Any other character string}{Evaluated as a \code{dplyr::mutate()}
#'       expression directly.}
#'   }
#' @param verbose Logical.  If \code{TRUE}, informational messages are printed.
#'   Default is \code{FALSE}.
#' @param memory_manage Integer memory-management flag controlling whether
#'   intermediate garbage collection is performed and whether certain columns are
#'   dropped to conserve RAM.  Default is \code{0L}.
#'
#' @return A data frame or tibble with the requested derived columns added.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @keywords internal
mutateFunction <- function(tbl, mutate_vars, verbose = FALSE, memory_manage = 0L) {
  # subfunctions

  ## case control

  mutateCaseControl <- function(tbl, memory_manage = 0L) {
    if ("casecontrol_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
      }
      tbl %>% # match
        dplyr::mutate(
          casecontrol_groupings = dplyr::case_when(
            CaseControl_k1 == 1 & Has_Dementia_k1 == 1 ~ 1L,
            CaseControl_k1 == 1 & Has_Dementia_k1 == 0 ~ 0L,
            TRUE ~ NA_integer_
          )
        )
    }
  }


  ## gender
  mutateGender <- function(tbl, memory_manage = 0L) {
    if ("gender_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            gender_groupings = dplyr::case_when(
              .data$male_k1 == 1 & .data$male_k2 == 1 ~ 2L,
              .data$male_k1 == 0 & .data$male_k2 == 0 ~ 0L,
              .data$male_k1 == 1 & .data$male_k2 == 0 ~ 1L,
              .data$male_k1 == 0 & .data$male_k2 == 1 ~ 1L,
              TRUE ~ NA_integer_
            )
          ) %>% dplyr::select(-"male_k2")
      } else {
        tbl %>% # match
          dplyr::mutate(
            gender_groupings = dplyr::case_when(
              .data$male_k1 == 1 & .data$male_k2 == 1 ~ 2L,
              .data$male_k1 == 0 & .data$male_k2 == 0 ~ 0L,
              .data$male_k1 == 1 & .data$male_k2 == 0 ~ 1L,
              .data$male_k1 == 0 & .data$male_k2 == 1 ~ 1L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  ## linkage type
  mutateLinkage <- function(tbl, memory_manage = 0L) {
    if ("linkagetype" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("patID_k1", "patID_k2", "matID_k1", "matID_k2")) %>%
          dplyr::mutate(linkagetype = dplyr::case_when(
            .data$same_matID == 1 & .data$same_patID == 1 ~ 11L,
            .data$same_matID == 0 & .data$same_patID == 0 ~ 00L,
            .data$same_matID == 1 & .data$same_patID == 0 ~ 10L,
            .data$same_matID == 0 & .data$same_patID == 1 ~ 01L,
            TRUE ~ NA_integer_
          )) %>%
          dplyr::select(-c("same_matID", "same_patID"))
      } else {
        tbl %>% # match
          dplyr::mutate(
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            linkagetype = dplyr::case_when(
              .data$same_matID == 1 & .data$same_patID == 1 ~ 11L,
              .data$same_matID == 0 & .data$same_patID == 0 ~ 00L,
              .data$same_matID == 1 & .data$same_patID == 0 ~ 10L,
              .data$same_matID == 0 & .data$same_patID == 1 ~ 01L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }

  mutateLinkage_same_patID <- function(tbl, memory_manage = 0L) {
    if ("same_patID" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            #  same_matID = dplyr::case_when(
            #    matID_k1 == matID_k2 ~ 1L,
            #    matID_k1 != matID_k2 ~ 0L,
            #    TRUE ~ NA_integer_
            #  ),
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          ) %>% dplyr::select(-c("patID_k1", "patID_k2", "matID_k1", "matID_k2"))
      } else {
        tbl %>% # match
          dplyr::mutate(
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  mutateLinkage_same_matID <- function(tbl, memory_manage = 0L) {
    if ("same_matID" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            #  same_matID = dplyr::case_when(
            #    matID_k1 == matID_k2 ~ 1L,
            #    matID_k1 != matID_k2 ~ 0L,
            #    TRUE ~ NA_integer_
            #  ),
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          ) %>% dplyr::select(-c("patID_k1", "patID_k2", "matID_k1", "matID_k2"))
      } else {
        tbl %>% # match
          dplyr::mutate(
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }

  mutateLinkage_any <- function(tbl, memory_manage = 0L) {
    if ("linkagetype" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("patID_k1", "patID_k2", "matID_k1", "matID_k2")) %>%
          dplyr::mutate(linkagetype = dplyr::case_when(
            .data$same_matID == 1 & .data$same_patID == 1 ~ 11L,
            .data$same_matID == 0 & .data$same_patID == 0 ~ 00L,
            .data$same_matID == 1 & .data$same_patID == 0 ~ 11L,
            .data$same_matID == 0 & .data$same_patID == 1 ~ 11L,
            TRUE ~ NA_integer_
          )) %>%
          dplyr::select(-c("same_matID", "same_patID"))
      } else {
        tbl %>% # match
          dplyr::mutate(
            same_matID = dplyr::case_when(
              .data$matID_k1 == .data$matID_k2 ~ 1L,
              .data$matID_k1 != .data$matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = dplyr::case_when(
              .data$patID_k1 == .data$patID_k2 ~ 1L,
              .data$patID_k1 != .data$patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            linkagetype = dplyr::case_when(
              .data$same_matID == 1 & .data$same_patID == 1 ~ 11L,
              .data$same_matID == 0 & .data$same_patID == 0 ~ 00L,
              .data$same_matID == 1 & .data$same_patID == 0 ~ 11L,
              .data$same_matID == 0 & .data$same_patID == 1 ~ 11L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  mutateCohort_19 <- function(tbl, memory_manage = 0L) {
    if ("cohort_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("BYr_k1", "BYr_k2")) %>%
          dplyr::mutate(
            cohort_groupings = dplyr::case_when(
              .data$cohort_k1 < .data$cohort_k2 ~ as.integer(.data$cohort_k1 * 100 + .data$cohort_k2),
              .data$cohort_k1 >= .data$cohort_k2 ~ as.integer(.data$cohort_k2 * 100 + .data$cohort_k1),
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("cohort_k1", "cohort_k2"))
      } else {
        tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              .data$cohort_k1 < .data$cohort_k2, paste(.data$cohort_k1, .data$cohort_k2, sep = "_"),
              paste(.data$cohort_k2, .data$cohort_k1, sep = "_")
            )
          )
      }
    }
  }
  mutateCohort_19flat <- function(tbl, memory_manage = 0L) {
    if ("cohort_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("BYr_k1", "BYr_k2")) %>%
          dplyr::mutate(
            cohort_groupings = dplyr::case_when(
              .data$cohort_k1 == 19 & .data$cohort_k2 == 19 ~ as.integer(1919),
              .data$cohort_k1 < 19 | .data$cohort_k2 < 19 ~ as.integer(1818),
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("cohort_k1", "cohort_k2"))
      } else {
        tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ), cohort_groupings = dplyr::case_when(
              .data$cohort_k1 == 19 & .data$cohort_k2 == 19 ~ as.integer(1919),
              .data$cohort_k1 < 19 | .data$cohort_k2 < 19 ~ as.integer(1818),
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  # cohort_century
  mutateCohortMatch <- function(tbl, memory_manage = 0L) {
    if ("cohort_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 >= 1600 & .data$BYr_k1 < 1700 ~ 16L,
              .data$BYr_k1 >= 1700 & .data$BYr_k1 < 1800 ~ 17L,
              .data$BYr_k1 >= 1800 & .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 & .data$BYr_k1 < 2000 ~ 19L,
              .data$BYr_k1 >= 2000 & .data$BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 >= 1600 & .data$BYr_k2 < 1700 ~ 16L,
              .data$BYr_k2 >= 1700 & .data$BYr_k2 < 1800 ~ 17L,
              .data$BYr_k2 >= 1800 & .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 & .data$BYr_k2 < 2000 ~ 19L,
              .data$BYr_k2 >= 2000 & .data$BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("BYr_k1", "BYr_k2")) %>%
          dplyr::mutate(
            cohort_groupings = dplyr::case_when(
              .data$cohort_k1 == .data$cohort_k2 ~ as.integer(1),
              .data$cohort_k1 != .data$cohort_k2 ~ as.integer(0),
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("cohort_k1", "cohort_k2"))
      } else {
        tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 >= 1600 & .data$BYr_k1 < 1700 ~ 16L,
              .data$BYr_k1 >= 1700 & .data$BYr_k1 < 1800 ~ 17L,
              .data$BYr_k1 >= 1800 & .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 & .data$BYr_k1 < 2000 ~ 19L,
              .data$BYr_k1 >= 2000 & .data$BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 >= 1600 & .data$BYr_k2 < 1700 ~ 16L,
              .data$BYr_k2 >= 1700 & .data$BYr_k2 < 1800 ~ 17L,
              .data$BYr_k2 >= 1800 & .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 & .data$BYr_k2 < 2000 ~ 19L,
              .data$BYr_k2 >= 2000 & .data$BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              .data$cohort_k1 == .data$cohort_k2, as.integer(1),
              as.integer(0)
            )
          )
      }
    }
  }

  # cohort_century
  mutateCohort <- function(tbl, memory_manage = 0L) {
    if ("cohort_groupings" %in% names(tbl)) {
      tbl # skip in outcome vare already present
    } else {
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 >= 1600 & .data$BYr_k1 < 1700 ~ 16L,
              .data$BYr_k1 >= 1700 & .data$BYr_k1 < 1800 ~ 17L,
              .data$BYr_k1 >= 1800 & .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 & .data$BYr_k1 < 2000 ~ 19L,
              .data$BYr_k1 >= 2000 & .data$BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 >= 1600 & .data$BYr_k2 < 1700 ~ 16L,
              .data$BYr_k2 >= 1700 & .data$BYr_k2 < 1800 ~ 17L,
              .data$BYr_k2 >= 1800 & .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 & .data$BYr_k2 < 2000 ~ 19L,
              .data$BYr_k2 >= 2000 & .data$BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("BYr_k1", "BYr_k2")) %>%
          dplyr::mutate(
            cohort_groupings = dplyr::case_when(
              .data$cohort_k1 < .data$cohort_k2 ~ as.integer(.data$cohort_k1 * 100 + .data$cohort_k2),
              .data$cohort_k1 >= .data$cohort_k2 ~ as.integer(.data$cohort_k2 * 100 + .data$cohort_k1),
              TRUE ~ NA_integer_
            )
          ) %>%
          dplyr::select(-c("cohort_k1", "cohort_k2"))
      } else {
        tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = dplyr::case_when(
              .data$BYr_k1 >= 1600 & .data$BYr_k1 < 1700 ~ 16L,
              .data$BYr_k1 >= 1700 & .data$BYr_k1 < 1800 ~ 17L,
              .data$BYr_k1 >= 1800 & .data$BYr_k1 < 1900 ~ 18L,
              .data$BYr_k1 >= 1900 & .data$BYr_k1 < 2000 ~ 19L,
              .data$BYr_k1 >= 2000 & .data$BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = dplyr::case_when(
              .data$BYr_k2 >= 1600 & .data$BYr_k2 < 1700 ~ 16L,
              .data$BYr_k2 >= 1700 & .data$BYr_k2 < 1800 ~ 17L,
              .data$BYr_k2 >= 1800 & .data$BYr_k2 < 1900 ~ 18L,
              .data$BYr_k2 >= 1900 & .data$BYr_k2 < 2000 ~ 19L,
              .data$BYr_k2 >= 2000 & .data$BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              .data$cohort_k1 < .data$cohort_k2, paste(.data$cohort_k1, .data$cohort_k2, sep = "_"),
              paste(.data$cohort_k2, .data$cohort_k1, sep = "_")
            )
          )
      }
    }
  }
  # Main
  if (memory_manage > 0L) {
    gc()
  }
  if (length(mutate_vars) == 0 || is.na(mutate_vars) || mutate_vars == 0) {
    tbl
    #   message(paste("mutateFunction 1"))
  } else if (mutate_vars == "gender_groupings_linkagetype" | mutate_vars == "gender_linkagetype") {
    tbl %>%
      mutateLinkage(memory_manage = memory_manage) %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mutate_vars == "linkagetype") {
    tbl %>% # match
      mutateLinkage(memory_manage = memory_manage)
  } else if (mutate_vars == "linkage_any") {
    tbl %>% # match
      mutateLinkage_any(memory_manage = memory_manage)
  } else if (mutate_vars == "same_patID") {
    tbl %>% # match
      mutateLinkage_same_patID(memory_manage = memory_manage)
  } else if (mutate_vars == "same_matID") {
    tbl %>% # match
      mutateLinkage_same_matID(memory_manage = memory_manage)
  } else if (mutate_vars == "gender_groupings") {
    tbl %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mutate_vars == "casecontrol_groupings") {
    tbl %>%
      mutateCaseControl(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_groupings_19") {
    tbl %>% mutateCohort_19(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_groupings_19flat") {
    tbl %>% mutateCohort_19flat(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_groupings_match") {
    tbl %>% mutateCohortMatch(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_groupings") {
    tbl %>% mutateCohort(memory_manage = memory_manage)
  } else if (mutate_vars == "gender_cohort_groupings_19") {
    tbl %>%
      mutateCohort_19(memory_manage = memory_manage) %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mutate_vars == "gender_cohort_groupings") {
    tbl %>%
      mutateCohort(memory_manage = memory_manage) %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_linkagetype_19") {
    tbl %>%
      mutateLinkage(memory_manage = memory_manage) %>%
      mutateCohort_19(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_linkagetype") {
    tbl %>%
      mutateLinkage(memory_manage = memory_manage) %>%
      mutateCohort(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_gender_linkagetype_19") {
    tbl %>%
      mutateLinkage(memory_manage = memory_manage) %>%
      mutateCohort_19(memory_manage = memory_manage) %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_gender_linkagetype") {
    tbl %>%
      mutateLinkage(memory_manage = memory_manage) %>%
      mutateCohort(memory_manage = memory_manage) %>%
      mutateGender(memory_manage = memory_manage)
  } else if (mode(mutate_vars) == "character") {
    mutate_call <- paste0("tbl %>% dplyr::mutate(", mutate_vars, ")")
    eval(parse(text = mutate_call))
  } else {
    tbl
  }
}


# creates the group bys
## I tried A TON OF THINGS, but the most predictably behaving version...
## was to create a large string via mapply
#' Summarize Outcomes Across Kinship Bins
#'
#' Builds and evaluates a \code{dplyr::summarize()} call that computes
#' within-bin statistics for each outcome variable using the specified
#' estimation function.  The summary always includes pair counts and
#' relatedness-bin metadata; per-outcome columns are appended by constructing
#' an R expression string via \code{mapply()} and evaluating it.
#'
#' @param tbl A (grouped or ungrouped) data frame containing the kinship-pair
#'   data for a single bin.
#' @param outcome_k1 Character vector of outcome variable names for the first
#'   member of each pair (the \code{_k1} side).  Must be the same length as
#'   \code{outcome_k2} and \code{outcome_functions}.
#' @param outcome_k2 Character vector of outcome variable names for the second
#'   member of each pair (the \code{_k2} side).
#' @param outcome_functions Character vector naming the statistical function to
#'   apply to each outcome pair.  Supported values include
#'   \code{"meanFunction"}, \code{"sdFunction"}, \code{"polychorFunction"},
#'   \code{"ml_polychorFunction"}, \code{"relriskFunction"},
#'   \code{"phi_both"}, \code{"phi_none"}, \code{"phi_one"},
#'   \code{"phi_k1_yes_k2_no"}, \code{"phi_k1_no_k2_yes"},
#'   \code{"phi_est"}, \code{"phi_se"}, \code{"phi_LL"}, \code{"phi_UL"},
#'   \code{"phi_ci"}, \code{"rr_exposed_cases"}, \code{"rr_unexposed_noncases"},
#'   \code{"rr_discordant"}, \code{"rr_exposed_noncases"},
#'   \code{"rr_unexposed_cases"}, \code{"rr_risk_exposed"},
#'   \code{"rr_risk_unexposed"}, \code{"rr_est"}, and any other function
#'   name that takes a single vector argument.
#' @param mitj Value of the mitochondrial-DNA relatedness flag for the current
#'   iteration (stored in the output as \code{mtdna}).
#' @param cnuk Value of the common-nuclear-environment flag for the current
#'   iteration (stored in the output as \code{cnu}).
#' @param range_maxi Upper bound of the additive-relatedness bin.
#' @param range_mini Lower bound of the additive-relatedness bin.
#' @param verbose Logical.  If \code{TRUE}, the constructed \code{summarize}
#'   call string is printed via \code{message()}.  Default is \code{FALSE}.
#' @param memory_manage Integer memory-management flag.  Values \code{> 1}
#'   suppress empirical \code{addRel} statistics to save RAM.  Default is
#'   \code{0L}.
#' @param skinny_summarize_call Logical.  If \code{TRUE} (default), the
#'   summary omits empirical \code{addRel} distribution statistics
#'   (\code{min}, \code{mean}, \code{median}, \code{max}) and the unique-N
#'   count.  Set to \code{FALSE} for a richer output.
#' @param SEN Logical.  If \code{TRUE}, an extended set of summary statistics
#'   (including unique individual counts and average dyads per ID) is computed.
#'   Default is \code{FALSE}.
#'
#' @return A one-row (or one-row-per-group) tibble containing the summary
#'   statistics for the current bin.
#'
#' @importFrom magrittr %>%
#' @keywords internal
summarizerFunction <- function(tbl,
                               outcome_k1,
                               outcome_k2,
                               outcome_functions,
                               mitj,
                               cnuk,
                               range_maxi,
                               range_mini,
                               verbose = FALSE, memory_manage = 0L,
                               skinny_summarize_call = TRUE,
                               SEN = FALSE) {
  if (length(outcome_k1) != length(outcome_functions)) {
    stop("The vectors of function names and variables must be the same length")
  }
  if (length(outcome_k1) != length(outcome_k2)) {
    stop("The vectors of outcome names must be the same length")
  }

  if (SEN == FALSE) {
    if (skinny_summarize_call == TRUE) {
      if (memory_manage > 1 || !("addRel" %in% names(tbl))) {
        summarize_call <- "tbl %>% dplyr::summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      } else {
        summarize_call <- "tbl %>% dplyr::summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      addRel_emp_min = try_NA(min(addRel,na.rm=TRUE)),
      addRel_emp_mean = try_NA(mean(addRel,na.rm=TRUE)),
      addRel_emp_median = try_NA(median(addRel,na.rm=TRUE)),
      addRel_emp_max = try_NA(max(addRel,na.rm=TRUE)),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      }
    } else {
      if (memory_manage > 1 || !("addRel" %in% names(tbl))) {
        summarize_call <- "tbl %>% dplyr::summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = dplyr::n_distinct(c(ID2,ID1)),
	    addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      } else {
        summarize_call <- "tbl %>% dplyr::summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = dplyr::n_distinct(c(ID2,ID1)),
	    addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      addRel_emp_min = try_NA(min(addRel,na.rm=TRUE)),
      addRel_emp_mean = try_NA(mean(addRel,na.rm=TRUE)),
      addRel_emp_median = try_NA(median(addRel,na.rm=TRUE)),
      addRel_emp_max = try_NA(max(addRel,na.rm=TRUE)),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      }
    }
  } else {
    summarize_call <- "tbl %>%
   dplyr::summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = dplyr::n_distinct(c(ID1, ID2)),
      avg_rows_per_id = n() / dplyr::n_distinct(c(ID1, ID2)),
      avg_dyads_per_id = (n() / (1 + doubleentered)) / dplyr::n_distinct(c(ID1, ID2)),
    #  icc_by_id = estimate_icc_latent_from_dyadic(outcome_var='USA_flag_10', method = 'mean'),
	    addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
  }


  summarize_parts <- mapply(
    function(var_k1,
             var_k2,
             fun) {
      if (fun == "polychorFunction") {
        paste0(
          var_k1, "_", fun, " = list(try_NA(polycor::polychor(",
          var_k1, "_k1,",
          var_k2, "_k2, std.err=TRUE)) %>%
                 {list(rho = try_NA(.$rho),
                       se = sqrt(try_NA(.$var)),
                       chisq = try_NA(.$chisq),
                       df = try_NA(.$df))})"
        )
      } else if (fun == "ml_polychorFunction") {
        paste0(var_k1, "_", fun, " = try_NA(polycor::polychor(", var_k1, "_k1,", var_k2, "_k2,ML=TRUE))")
      } else if (fun == "relriskFunction") {
        paste0(
          var_k1, "_", fun, " = list(try_NA(relriskFunction(",
          var_k1, "_k1, ",
          var_k2, "_k2)) %>%
             {list(rr = try_NA(.[1]),
                   LL = try_NA(.[2]),
                   UL = try_NA(.[3]))})"
        )
      } else if (fun %in% c("phi_both", "rr_exposed_cases", "rr_a")) {
        paste0(var_k1, "_", fun, "= try_NA(sum(", var_k1, "_k1 == 1 & ", var_k2, "_k2 == 1, na.rm=TRUE))")
      } else if (fun %in% c("phi_none", "rr_unexposed_noncases", "rr_d")) {
        paste0(var_k1, "_", fun, "= try_NA(sum(", var_k1, "_k1 == 0 & ", var_k2, "_k2 == 0, na.rm=TRUE))")
      } else if (fun %in% c("phi_one", "rr_discordant", "rr_b_plus_c")) {
        paste0(var_k1, "_", fun, "= try_NA(sum((", var_k1, "_k1 == 1 & ", var_k2, "_k2 == 0) | (", var_k1, "_k1 == 0 & ", var_k2, "_k2 == 1), na.rm = TRUE))")
      } else if (fun %in% c("phi_k1_yes_k2_no", "rr_exposed_noncases", "rr_b")) {
        paste0(var_k1, "_", fun, " = try_NA(sum(", var_k1, "_k1 == 1 & ", var_k2, "_k2 == 0, na.rm = TRUE))")
      } else if (fun %in% c("phi_k1_no_k2_yes", "rr_unexposed_cases", "rr_c")) {
        paste0(var_k1, "_", fun, " = try_NA(sum(", var_k1, "_k1 == 0 & ", var_k2, "_k2 == 1, na.rm = TRUE))")
      } else if (fun == "phi_est") {
        paste0(
          var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
          var_k1, "_phi_both, f01= ",
          var_k1, "_phi_one, f10=",
          var_k1, "_phi_one, f00 = ",
          var_k1, "_phi_none))[1]"
        )
      } else if (fun == "phi_se") {
        paste0(
          var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
          var_k1, "_phi_both, f01= ",
          var_k1, "_phi_one,f10=",
          var_k1, "_phi_one,f00 = ",
          var_k1, "_phi_none))[2]"
        )
      } else if (fun == "phi_LL") {
        paste0(
          var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
          var_k1, "_phi_both, f01= ",
          var_k1, "_phi_one,f10=",
          var_k1, "_phi_one,f00 = ",
          var_k1, "_phi_none))[3]"
        )
      } else if (fun == "phi_UL") {
        paste0(
          var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
          var_k1, "_phi_both, f01= ",
          var_k1, "_phi_one,f10=",
          var_k1, "_phi_one,f00 = ",
          var_k1, "_phi_none))[4]"
        )
      } else if (fun == "phi_ci") {
        paste0(
          var_k1, "_", fun, " = list(try_NA(ci.phi(.05,f11= ",
          var_k1, "_phi_both, f01= ",
          var_k1, "_phi_one,f10=",
          var_k1, "_phi_one,f00 = ",
          var_k1, "_phi_none))  %>%
             {list(phi = try_NA(.$Estimate),
                   se = try_NA(.$SE),
                   LL = try_NA(.$LL),
                   UL = try_NA(.$UL))})"
        )
      } else if (fun == "rr_risk_exposed") {
        paste0(
          var_k1, "_", fun_k1, " = try_NA(ifelse((",
          var_k1, "_rr_a + ", var_k1, "_rr_b) > 0, ",
          var_k1, "_rr_a / (", var_k1, "_rr_a + ", var_k1, "_rr_b), ",
          "NA_real_))"
        )
      } else if (fun == "rr_risk_unexposed") {
        paste0(
          var_k1, "_", fun, " = try_NA(ifelse((",
          var_k1, "_rr_c + ", var_k1, "_rr_d) > 0, ",
          var_k1, "_rr_c / (", var_k1, "_rr_c + ", var_k1, "_rr_d), ",
          "NA_real_))"
        )
      } else if (fun == "rr_est") {
        paste0(
          var_k1, "_", fun, " = try_NA(ifelse(",
          var_k1, "_rr_risk_unexposed > 0, ",
          var_k1, "_rr_risk_exposed / ", var, "_rr_risk_unexposed, ",
          "NA_real_))"
        )
      } else {
        paste0(var_k1, "_", fun, " = try_NA(", fun, "(", var_k1, "))")
      }
    }, outcome_k1, outcome_k2,
    outcome_functions,
    SIMPLIFY = FALSE
  )
  # if unnesting variable is needed
  unnest_parts <- c()

  if ("polychorFunction" %in% outcome_functions) {
    # get vars to unlist
    # doesn't like when there are multiple polychors
    var_polychor_unnest <- outcome_k1[outcome_functions == "polychorFunction"]
    unnest_parts <- c(
      unnest_parts,
      mapply(function(var) {
        paste0("unnest_wider(", var, "_polychorFunction, names_sep = '_')")
      }, var_polychor_unnest, SIMPLIFY = FALSE)
    )
  }

  if ("phi_ci" %in% outcome_functions) {
    var_phi_unnest <- outcome_k1[outcome_functions == "phi_ci"]
    unnest_parts <- c(
      unnest_parts,
      mapply(function(var) {
        paste0("unnest_wider(", var, "_phi_ci, names_sep = '_')")
      }, var_phi_unnest, SIMPLIFY = FALSE)
    )
  }
  if ("relriskFunction" %in% outcome_functions) {
    var_rr_unnest <- outcome_k1[outcome_functions == "relriskFunction"]
    unnest_parts <- c(
      unnest_parts,
      mapply(function(var) {
        paste0("unnest_wider(", var, "_relriskFunction, names_sep = '_')")
      }, var_rr_unnest, SIMPLIFY = FALSE)
    )
  }

  # combine all the parts into a single call string
  if (length(unnest_parts) > 0) {
    unnest_call <- paste0(" %>% ", paste(unnest_parts, collapse = " %>% "))
  } else {
    unnest_call <- ""
  }


  summarize_call <- paste0(
    summarize_call,
    paste(summarize_parts, collapse = ", "), ")", unnest_call
  )
  if (verbose) {
    message(summarize_call)
  }

  # evaluate the constructed function call and return result
  eval(parse(text = summarize_call))
}


#' Compute the Mean, Ignoring Missing Values
#'
#' A thin wrapper around \code{base::mean()} with \code{na.rm = TRUE}.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value: the arithmetic mean of the non-missing
#'   elements of \code{x}.
#'
#' @keywords internal
meanFunction <- function(x) {
  mean(x, na.rm = TRUE)
}


#' Compute the Standard Deviation, Ignoring Missing Values
#'
#' A thin wrapper around \code{stats::sd()} with \code{na.rm = TRUE}.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value: the sample standard deviation of the
#'   non-missing elements of \code{x}.
#'
#' @importFrom stats sd
#' @keywords internal
sdFunction <- function(x) {
  sd(x, na.rm = TRUE)
}

#' Compute the 25th Percentile, Ignoring Missing Values
#'
#' A thin wrapper around \code{stats::quantile()} that returns the 25th
#' percentile (first quartile) with \code{na.rm = TRUE}.
#'
#' @param x A numeric vector.
#' @return A single numeric value: the 25th percentile of \code{x}.
#' @importFrom stats quantile
#' @keywords internal

q25Function <- function(x) {
  quantile(x, na.rm = TRUE, probs = .25)
}

#' Compute the 75th Percentile, Ignoring Missing Values
#'
#' A thin wrapper around \code{stats::quantile()} that returns the 75th
#' percentile (third quartile) with \code{na.rm = TRUE}.
#'
#' @param x A numeric vector.
#' @return A single numeric value: the 75th percentile of \code{x}.
#' @importFrom stats quantile
#' @keywords internal

q75Function <- function(x) {
  quantile(x, na.rm = TRUE, probs = .75)
}

#' Compute the Median (50th Percentile), Ignoring Missing Values
#'
#' A thin wrapper around \code{stats::quantile()} that returns the 50th
#' percentile (median) with \code{na.rm = TRUE}.
#'
#' @param x A numeric vector.
#' @return A single numeric value: the median of \code{x}.
#' @importFrom stats quantile
#' @keywords internal

q50Function <- function(x) {
  quantile(x, na.rm = TRUE, probs = .5)
}

#' Compute the Relative Risk for a Binary Outcome Pair
#'
#' Constructs a 2 \eqn{\times} 2 contingency table from two binary vectors
#' and estimates the relative risk (and confidence interval) using
#' \code{DescTools::RelRisk()}.
#'
#' @param k1 A binary numeric (or integer) vector for the first kin member
#'   (values expected to be 0 or 1).
#' @param k2 A binary numeric (or integer) vector for the second kin member,
#'   the same length as \code{k1}.
#' @param conf.level Numeric confidence level for the interval.  Default is
#'   \code{0.95}.
#' @param method Character string passed to \code{DescTools::RelRisk()}.
#'   Default is \code{"score"}.
#'
#' @return A named numeric vector with the relative risk estimate, lower
#'   confidence bound, and upper confidence bound (as returned by
#'   \code{DescTools::RelRisk()}).
#'
#' @keywords internal

relriskFunction <- function(k1, k2, conf.level = 0.95, method = "score") {
  rr_tbl <- table(
    factor(k1, levels = c(1, 0)),
    factor(k2, levels = c(1, 0))
  )
  DescTools::RelRisk(
    rr_tbl,
    conf.level = conf.level,
    method = method
  )
}


#' Build an Input File Path for a Kinship Bin
#'
#' Assembles the expected file path for a CSV file that stores relatedness-bin
#' data, following the naming convention used throughout this package.
#'
#' @param data_path Character string giving the root data directory (may be an
#'   empty string \code{""} if the path is relative to the working directory).
#' @param df_foldername Character string: name of the dataset/folder.
#' @param binwidth_cha Character string representation of the bin width (e.g.
#'   \code{"10"} for a 10 \% bin).
#' @param mit Value of the mitochondrial-DNA relatedness flag (\code{0} or
#'   \code{1}).
#' @param range_min Lower bound of the additive-relatedness range for this bin.
#' @param range_max Upper bound of the additive-relatedness range for this bin.
#'
#' @return A single character string giving the full path to the expected input
#'   CSV file.
#'
#' @keywords internal
make_input_file <- function(data_path,
                            df_foldername,
                            binwidth_cha, mit,
                            range_min,
                            range_max) {
  paste0(data_path, "data/", df_foldername, "_", binwidth_cha, "/df_mt", mit, "_r", range_min, "-r", range_max, ".csv")
}


#' Read a Kinship Bin CSV File
#'
#' Reads a single kinship-bin CSV file using \code{data.table::fread()},
#' optionally dropping specified columns, and returns the result as a tibble.
#' Errors from \code{fread()} are caught and \code{NA} is returned instead.
#'
#' @param input_file Character string: full path to the CSV file to read.
#' @param drop_variables Character vector of column names to drop when reading
#'   the file.  Default is \code{c("mitRel")}.
#' @param verbose Logical.  If \code{TRUE}, a message reporting the number of
#'   rows read is emitted.  Default is \code{FALSE}.
#'
#' @return A data frame (tibble-compatible) with the contents of
#'   \code{input_file}, or \code{NA} if reading fails.
#'
#' @importFrom magrittr %>%
#' @keywords internal
read_kinbin <- function(input_file, drop_variables = c("mitRel"), verbose = FALSE) {
  dataRelatedPair_merge <- try_NA(data.table::fread(input_file,
    header = TRUE,
    #  drop vars to slim
    drop = drop_variables
  )) %>%
    #   dplyr::mutate(addRel = round(addRel, digits = 4)) %>% # can be dropped?
    suppressWarnings()

  if (verbose == TRUE) {
    message(paste0(input_file, "had ", nrow(dataRelatedPair_merge), " rows"))
  }

  return(dataRelatedPair_merge)
}
