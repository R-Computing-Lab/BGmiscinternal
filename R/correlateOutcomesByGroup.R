#------------------------------------------------------------------------------
# Author Mason Garrison
# Date: 815
# Filename: functions_CorrelateOutcomesByGroup
# Purpose: this code calcuates the correlation between the outcomes for each kin group, groups by mtdna, cnu, and bins of R

#------------------------------------------------------------------------------

## First are the helper functions
options(scipen = 10, digits = 11)
## passes subset if slice_1000 is true, otherwise passes entire thing



# here's the mega function
#' Correlate Outcomes by Kinship Group
#'
#' The main workhorse function of this package.  Iterates over bins of
#' additive-genetic relatedness (R), mitochondrial-DNA (mtDNA) relatedness,
#' and common nuclear-environment (CNU) status; reads per-bin data files;
#' optionally double-enters each dyadic data set; applies user-specified
#' mutation steps; groups the data; and computes correlation/association
#' statistics for each specified outcome pair.  Results are written
#' incrementally to a CSV file.
#'
#' @param df_foldername Character string.  Name of the data folder/dataset.
#'   Default is \code{"longevity_skinny_matpat"}.
#' @param binwidth_num Numeric vector of bin half-widths (e.g. \code{0.10} for
#'   a \eqn{\pm 10\%} window around each relatedness center).  Must be the
#'   same length as \code{binwidth_cha}.  Default is \code{c(0.1, 0.05)}.
#' @param binwidth_cha Character vector of bin-width labels used in file names
#'   (e.g. \code{"10"} corresponds to \code{0.10}).  Default is
#'   \code{c("10", "05")}.
#' @param kin_degree_max Integer.  Maximum kinship degree to process
#'   (\eqn{2^{-\text{degree}}} gives the relatedness center).
#'   Default is \code{12}.
#' @param kin_degree_min Integer.  Minimum kinship degree.  Default is
#'   \code{0}.
#' @param max_kin_per_bin Numeric.  Bins containing more rows than this value
#'   are skipped.  Default is \code{8.7e9} (effectively no limit).
#' @param cnu Integer vector of CNU values to loop over (typically
#'   \code{c(1, 0)}).  Default is \code{c(1, 0)}.
#' @param mit Integer vector of mtDNA relatedness values to loop over
#'   (typically \code{c(1, 0)}).  Default is \code{c(1, 0)}.
#' @param doubleentered Logical.  If \code{TRUE} (default) the dyadic data are
#'   double-entered (each pair appears twice with \emph{k1}/\emph{k2} roles
#'   swapped) to ensure symmetry.
#' @param slice_1000 Logical.  If \code{TRUE}, only the first 1 000 rows of
#'   each bin are processed (useful for testing).  Default is \code{FALSE}.
#' @param cleanup Logical.  If \code{TRUE} (default), temporary \code{.RDS}
#'   scratch files created during memory-managed runs are deleted on exit.
#' @param drop_variables Character vector of column names to drop when reading
#'   each bin file.  Default drops several quantile columns.
#' @param outcome_vars Character vector of outcome variable names (base names
#'   without the \code{_k1}/\code{_k2} suffix).  When supplied,
#'   \code{outcome_k1} and \code{outcome_k2} are both set to this vector.
#'   At least one of \code{outcome_vars}, \code{outcome_k1}, or
#'   \code{outcome_k2} must be non-\code{NULL}.  Default is \code{NULL}.
#' @param outcome_k1 Character vector of \emph{k1}-side outcome column names.
#'   Ignored when \code{outcome_vars} is non-\code{NULL}.  Default is
#'   \code{NULL}.
#' @param outcome_k2 Character vector of \emph{k2}-side outcome column names.
#'   Same length as \code{outcome_k1}.  Ignored when \code{outcome_vars} is
#'   non-\code{NULL}.  Default is \code{NULL}.
#' @param outcome_functions Character vector of function names to apply to each
#'   outcome pair.  See \code{\link{summarizerFunction}} for supported values.
#'   Default is \code{c("meanFunction", "meanFunction", "meanFunction",
#'   "meanFunction", "polychorFunction", "ml_polychorFunction",
#'   "polychorFunction", "ml_polychorFunction")}.
#' @param grouping_vars Character vector of additional grouping variable names,
#'   or \code{NA} for no additional grouping.  Default is \code{NA}.
#' @param grouping_filename Character string used in the output file name to
#'   identify the grouping scheme, or \code{NULL} for none.  Default is
#'   \code{NULL}.
#' @param verbose Logical.  If \code{TRUE} (default), progress messages are
#'   emitted at each major loop iteration.
#' @param mutate_vars Character string (or \code{NULL}) naming the mutation
#'   scheme to apply before summarising.  Passed to
#'   \code{\link{mutateFunction}}.  Default is \code{NULL}.
#' @param memory_manage Integer flag controlling memory-optimization level:
#'   \code{0L} = no optimization; \code{1L} = write/read RDS scratch files;
#'   \code{2L} = additionally use \pkg{data.table}/\pkg{tidyft} operations.
#'   Default is \code{0L}.
#' @param data_path Character string giving the root path to the data
#'   directory.  Default is \code{""} (current working directory).
#' @param file_path_stem Character string giving the stem of the output CSV
#'   file path.  Default is \code{"U:/IRB_00143000/mtdna/aim1_cor/cor_"}.
#' @param skinny_summarize_call Logical.  If \code{TRUE} (default), the
#'   summary output omits empirical \code{addRel} distribution statistics.
#' @param age_filter Logical.  If \code{TRUE}, rows where both kin members are
#'   younger than \code{min_age} are removed before analysis.  Default is
#'   \code{FALSE}.
#' @param min_age Numeric.  Minimum age threshold used when
#'   \code{age_filter = TRUE}.  Default is \code{0}.
#' @param SEN Logical.  If \code{TRUE}, the summariser computes additional
#'   columns (unique individual counts, average dyads per ID).  Default is
#'   \code{FALSE}.
#'
#' @return Called primarily for its side effect of writing a CSV file to
#'   \code{file_path_stem}.  Returns \code{NULL} invisibly.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom dplyr filter across select
#' @importFrom tidyselect where
#' @importFrom data.table rbindlist fwrite
#' @importFrom readr write_rds read_rds
#' @importFrom utils read.csv
#' @export
correlateOutcomesByGroup <- function(df_foldername = "longevity_skinny_matpat",
                                     binwidth_num = c(.1, .05),
                                     binwidth_cha = c("10", "05"),
                                     kin_degree_max = 12,
                                     kin_degree_min = 0,
                                     max_kin_per_bin = 8.7 * 10^9, # 10^7 is 10 million
                                     cnu = c(1, 0),
                                     mit = c(1, 0),
                                     doubleentered = TRUE,
                                     slice_1000 = FALSE,
                                     cleanup = TRUE,
                                     drop_variables = c(
                                       "mitRel", "USA_quantile_k1", "USA_quantile_k2",
                                       "SWE_quantile_k1", "SWE_quantile_k2"
                                     ),
                                     # descriptive_vars=c("male","age"),
                                     #  outcome_vars=c("USA_flag_10",
                                     #                "USA_flag_15"),
                                     # alternative is var name and function to call
                                     # collumn of functions to call
                                     outcome_vars = NULL, #c(
                                     #   "male", "age",
                                     #   "USA_flag_10",
                                     #   "USA_flag_15",
                                     #    "USA_flag_10",
                                     #    "USA_flag_10",
                                     #    "USA_flag_15",
                                     #    "USA_flag_15"
                                     #  ),
                                     outcome_k1 = NULL,
                                     outcome_k2 = NULL,
                                     outcome_functions = c(
                                       "meanFunction",
                                       "meanFunction",
                                       "meanFunction",
                                       "meanFunction",
                                       "polychorFunction",
                                       "ml_polychorFunction",
                                       "polychorFunction",
                                       "ml_polychorFunction"
                                     ),
                                     grouping_vars = NA,
                                     grouping_filename = NULL,
                                     verbose = TRUE,
                                     mutate_vars = NULL,
                                     memory_manage = 0L,
                                     data_path = "",
                                     file_path_stem = "U:/IRB_00143000/mtdna/aim1_cor/cor_",
                                     skinny_summarize_call = TRUE,
                                     age_filter = FALSE, min_age = 0,
                                     SEN=FALSE) {
  # potential options
  # expand.grid on unique values and filter via a loop
  ## make sure to use verbose to give folks a sense of what they;re asking for
  ## include a check to estimate time?
  # groupby
  options(scipen = 10, digits = 11)
  if (memory_manage > 0L) {
    gc()
  }

  if(is.null(outcome_vars) && is.null(outcome_k1) && is.null(outcome_k2)){
    stop("At least one of the following inputs need to be provided: outcome_vars, outcome_k1, and outcome_k2")
  }

  if(!is.null(outcome_vars)){
    # overwriting
    outcome_k2 <-  outcome_k1 <- outcome_vars
  }

  if (verbose) {
    message("starting checks")
  }
  # Checks
  ## Is the bin length the same?
  if (length(binwidth_num) != length(binwidth_cha)) {
    # return error
    stop(paste("The length of binwidth_num and binwidth_cha don't match.
         binwidth_num has ", length(binwidth_num), "but binwidth_cha has ", length(binwidth_cha)))
  }
  ## Is range of degrees viable?
  if (!(mode(kin_degree_max) %in% c("numeric", "integer"))) {
    stop("kin_degree_max isn't a number")
  }
  if (!(mode(kin_degree_min) %in% c("numeric", "integer"))) {
    stop("kin_degree_min isn't a number")
  }
  if (kin_degree_max < 1 | kin_degree_max < kin_degree_min | kin_degree_max < 0) {
    stop("kin_degree_min or max isn't workable value. Either the value is too small or min is larger than max")
  }
  if (verbose) {
    message("ending checks")
  }
  # loop by bin width
  if (verbose) {
    message("bin width")
  }
  for (q in 1:length(binwidth_num)) {
    # select the bin for the loop
    bin_width_num_q <- binwidth_num[q]
    bin_width_cha_q <- binwidth_cha[q]

    if (verbose) {
      message(paste0("bin width ", bin_width_num_q))
    }
    # create relatedness center values
    ## these don't diff by bin, but it made more sense to keep them nearby
    addRel_center <- 2^(0:(-kin_degree_max)) # bin center

    addRel_maxs_temp <- addRel_center * (1 + bin_width_num_q)
    # inclusive
    addRel_mins_temp <- addRel_center * (1 - bin_width_num_q) # min_bin

    # this is supposed to have one of each
    addRel_real_maxs <- addRel_mins_temp[-length(addRel_mins_temp)]
    addRel_real_mins <- addRel_maxs_temp[-1]

    addRel_maxs <- c(
      1.5, sort(c(addRel_real_maxs, addRel_maxs_temp),
                decreasing = TRUE
      ),
      addRel_mins_temp[length(addRel_mins_temp)]
    )
    addRel_mins <- c(
      addRel_maxs_temp[1],
      sort(c(addRel_real_mins, addRel_mins_temp),
           decreasing = TRUE
      ), 0
    )

    kin_degree <- c(0, 0:(kin_degree_max + 1))

    # Make csv name
    ## is there a grouping variable
    if (is.null(grouping_filename) || length(grouping_filename) == 0 || is.na(grouping_filename)) {
      file_path_txt <- paste0(file_path_stem, df_foldername, "_", binwidth_cha[q], ".csv")
    } else {
      file_path_txt <- paste0(file_path_stem, grouping_filename, "_", df_foldername, "_", binwidth_cha[q], ".csv")
    }

    if (verbose) {
      message(file_path_txt)
    }

    if (verbose) {
      message("starting genetic loop")
    }
    # loop by genetic relatedness
    for (i in length(addRel_maxs):1) {
      range_max <- addRel_maxs[i]
      range_min <- addRel_mins[i]
      if (verbose) {
        message(range_min)
      }
      # loop by mit
      if (verbose) {
        message("starting mit loop")
      }
      for (j in 1:length(mit)) {
        # craft the file name

        input_file <- make_input_file(
          data_path = data_path,
          df_foldername = df_foldername,
          binwidth_cha = binwidth_cha[q],
          mit = mit[j],
          range_min = range_min,
          range_max = range_max
        )


        if (verbose) {
          message(input_file)
        }

        # check if file exists
        ## if not exist, skip
        if (!file.exists(input_file)) {
          # message("Your dumb file doesn't exist human")
          message(paste("Missing input file:", input_file))
          dataRelatedPair_merge <- NULL
          #     next
        } else {
          dataRelatedPair_merge <- read_kinbin(
            input_file = input_file,
            verbose = verbose,
            drop_variables = drop_variables
          )

          # add filter
          if (age_filter == TRUE) {
            full_rows <- nrow(dataRelatedPair_merge)

            if (all(c("age_k1", "age_k2") %in% names(dataRelatedPair_merge))) {
              dataRelatedPair_merge <- dataRelatedPair_merge %>%
                tidyr::drop_na("age_k1", "age_k2") %>%
                dplyr::filter(.data$age_k1 >= min_age & .data$age_k2 >= min_age)
            } else if ("age" %in% names(dataRelatedPair_merge)) {
              dataRelatedPair_merge <- dataRelatedPair_merge %>%
                tidyr::drop_na("age") %>%
                dplyr::filter(.data$age >= min_age)
            }
            if (verbose == TRUE) {
              message(paste(
                i, addRel_mins[i],
                "Filtering people who didn't live to",
                min_age, " ", full_rows - nrow(dataRelatedPair_merge), "removed"
              ))
            }
          }

          if (memory_manage > 0L) {
            #    gc()
          }

          if (memory_manage > 1L) {
            # round to 2 digits
            #  if(is.data.table(dataRelatedPair_merge)){
            #   dataRelatedPair_merge <- dataRelatedPair_merge %>% tidyft::mutate_vars(is.numeric,round, 2)
            #   } else {
            dataRelatedPair_merge <- dataRelatedPair_merge %>%
              dplyr::mutate(across(where(is.numeric), round, 2)) # %>%
            # convert_to_integer(memory_manage =  memory_manage)
            #    dataRelatedPair_merge <- convert_to_integer(dataRelatedPair_merge,memory_manage =  memory_manage)

            #    print(head(dataRelatedPair_merge))
            #  }
            gc()
          } else if (memory_manage > 0L) {
            # if(is.data.table(dataRelatedPair_merge)){
            #   dataRelatedPair_merge <- dataRelatedPair_merge %>% tidyft::mutate_vars(is.numeric,round, 3)
            #  } else {
            # possible optimization
            dataRelatedPair_merge <- dataRelatedPair_merge %>%
              dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), round, 3))

            # %>%
            # dataRelatedPair_merge <- convert_to_integer(dataRelatedPair_merge,memory_manage =  memory_manage)

            #  }

            gc()
          }
        }

        if (length(dataRelatedPair_merge) == 0 || nrow(dataRelatedPair_merge) > max_kin_per_bin || (length(dataRelatedPair_merge) == 1 && is.na(dataRelatedPair_merge))) {
          # this loop skips the file if it is too big, otherwise the loop fails
          # and gives unhappy warning and fails without cleaning up
          message(paste(
            i, addRel_mins[i],
            "was skipped because it was ",
            if(length(dataRelatedPair_merge) == 0){"empty"}else{paste0(nrow(dataRelatedPair_merge), "which is bigger than ", max_kin_per_bin)}
          ))
          remove(dataRelatedPair_merge)
          if (memory_manage > 0L) {
            gc()
          }
          next
        } else {
          if (doubleentered == TRUE) {
            # double entered
            if (verbose == TRUE) {
              message("Double entering the data")
            } # double entering is the obvious change to make to optimize
            if ("ID1" %in% names(dataRelatedPair_merge) && "ID2" %in% names(dataRelatedPair_merge)) {
              dxlist <- c(
                "ID1", "ID2", # intentional ordering
                "addRel", # "mitRel",
                "cnuRel",
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k2")],
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k1")]
              )

              dxlist_main <- c(
                "ID2", "ID1", # intentional ordering
                "addRel", # "mitRel",
                "cnuRel",
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k1")],
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k2")]
              )
            } else { # if IDs not there, this is a painlessway to reduce the size of the matrix
              dxlist <- c( # "ID1","ID2", #intentional ordering
                "addRel", # "mitRel",
                "cnuRel",
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k2")],
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k1")]
              )
              dxlist_main <- c(
                # "ID2", "ID1", # intentional ordering
                "addRel", # "mitRel",
                "cnuRel",
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k1")],
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k2")]
              )
            }
            if (verbose) {
              print(dxlist)
            }
            # remove addRel
            if ("addRel" %notin% names(dataRelatedPair_merge)) {
              dxlist <- dxlist[!dxlist %in% "addRel"]
              dxlist_main <- dxlist_main[!dxlist_main %in% "addRel"]
            }

            if ("cnuRel" %notin% names(dataRelatedPair_merge)) {
              dxlist <- dxlist[!dxlist %in% "cnuRel"]
              dxlist_main <- dxlist_main[!dxlist_main %in% "cnuRel"]
            }
            if(verbose){
              print(names(dataRelatedPair_merge[, .SD, .SDcols = dxlist_main]))
              print(names(dataRelatedPair_merge[, .SD, .SDcols = dxlist]))
            }
            dataRelatedPair_merge <- data.table::rbindlist(
              list(
                dataRelatedPair_merge[, .SD, .SDcols = dxlist_main],
                dataRelatedPair_merge[, .SD, .SDcols = dxlist]
              ),
              use.names = FALSE
            )
            gc()
            remove(dxlist)
            remove(dxlist_main)

          }
        }
        # loop by common environment
        if (verbose == TRUE) {
          message("loop by common environment")
        }
        for (k in 1:length(cnu)) {
          current_nrows <- nrow(dataRelatedPair_merge %>% dplyr::filter(.data$cnuRel == cnu[k]))

          if (verbose) {
            message(paste("number of rows for", cnu[k], current_nrows))
          }
          if (current_nrows == 0) {
            # skip if no cnu
            if (memory_manage > 0L) {
              gc()
            }
            next
          } else {
            if (verbose) {
              #   message(head(dataRelatedPair_merge))
            }
            #  message("here")
            # sometimes the nested functions can't find the looping index
            cnuk <- cnu[k]
            mitj <- mit[j]

            if (memory_manage > 1L) {
              gc()
              require(tidyft) # nolint: undesirable_function_linter
              # this writes each iteration to disk the idea being that you loop thru all the groups
              ## not great solution but... it does let you use mutate once

              dataRelatedPair_merge %>%
                mutateFunction(
                  mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage
                ) %>%
                convert_to_integer(memory_manage = memory_manage) %>%
                readr::write_rds("dataRelatedPair_merge.RDS")

              dataRelatedPair_merge <- tidyft::setDT(dataRelatedPair_merge) %>% dplyr::filter(.data$cnuRel == cnuk)
              gc()
              if (length(grouping_vars) == 1 && !is.na(grouping_vars) && grouping_vars %in% names(dataRelatedPair_merge)) {
                grouping_loop <- unique(dataRelatedPair_merge[[grouping_vars]])
              } else {
                grouping_loop <- 1
              }
              for (g in length(grouping_loop)) {
                if (length(grouping_loop) == 1 && grouping_loop == 1) {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- readr::read_rds("dataRelatedPair_merge.RDS") %>%
                    tidyft::setDT() %>%
                    dplyr::filter(.data$cnuRel == cnuk)
                  gc()
                } else {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- readr::read_rds("dataRelatedPair_merge.RDS") %>%
                    tidyft::setDT(dataRelatedPair_merge) %>%
                    dplyr::filter(.data[[grouping_vars]] == grouping_loop[g], .data$cnuRel == cnuk)
                  gc()
                }

                temp <- dataRelatedPair_merge %>%
                  sliceFunction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                  #    mutateFunction(mutate_vars=mutate_vars,verbose=verbose, memory_manage = memory_manage)%>%
                  group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  summarizerFunction(
                    outcome_k1,
                    outcome_k2,
                    outcome_functions,
                    cnuk = cnuk, mitj = mitj,
                    range_maxi = range_max,
                    range_mini = range_min, verbose = verbose,
                    memory_manage = memory_manage,
                    skinny_summarize_call = skinny_summarize_call,
                    SEN=SEN
                  )
                gc()
              }
              gc()
              # is needed because we might loop across multiple cnus
              dataRelatedPair_merge <- readr::read_rds("dataRelatedPair_merge.RDS")
            } else if (memory_manage > 0L) { # clean up the this to be one variable with 3 levels, character is fine, 0L
              gc()

              # this writes each iteration to disk the idea being that you loop thru all the groups
              ## not great solution but... it does let you use mutate once
              # print("353")
              dataRelatedPair_merge %>%
                # mutateFunction(
                #    mutate_vars = mutate_vars, verbose = verbose,  memory_manage= memory_manage
                # ) %>% convert_to_integer(memory_manage =  memory_manage) %>%
                readr::write_rds("dataRelatedPair_merge.RDS")

              dataRelatedPair_merge %>%
                mutateFunction(
                  mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage
                ) %>%
                dplyr::filter(.data$cnuRel == cnuk) %>%
                convert_to_integer(memory_manage = memory_manage) %>%
                dplyr::select(-"cnuRel") %>%
                readr::write_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS"))


              dataRelatedPair_merge <- dataRelatedPair_merge %>% dplyr::filter(.data$cnuRel == cnuk)
              gc()
              # only one grouping variable
              if (length(grouping_vars) == 1 && !is.na(grouping_vars) && grouping_vars %in% names(dataRelatedPair_merge)) {
                grouping_loop <- unique(dataRelatedPair_merge[[grouping_vars]])
              } else if (FALSE && length(grouping_vars) == 2 && grouping_vars %in% names(dataRelatedPair_merge)) {
                # obviously this needs to be generalized
                grouping_loop_grid <- expand.grid(
                  unique(dataRelatedPair_merge[[grouping_vars[1]]]),
                  unique(dataRelatedPair_merge[[grouping_vars[2]]])
                )
                # duplicated, check it out
                grouping_loop <- nrow(grouping_loop_grid)
              } else {
                grouping_loop <- 1
              }
              for (g in 1:length(grouping_loop)) {
                if (length(grouping_loop) == 1 && grouping_loop == 1) {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- readr::read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) # %>% # this can be managed
                  #  dplyr::filter(cnuRel == cnuk)
                  gc()
                } else if (FALSE) { # &exists("grouping_loop_grid")){
                  # if there's a grid
                  ## not complete
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- readr::read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) %>%
                    dplyr::filter(.data[[grouping_vars[1]]] == grouping_loop_grid[1, g], .data[[grouping_vars[2]]] == grouping_loop_grid[2, g])
                  gc()
                } else {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- readr::read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) %>%
                    dplyr::filter(.data[[grouping_vars]] == grouping_loop[g])
                  gc()
                }
                temp <- dataRelatedPair_merge %>%
                  sliceFunction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                  mutateFunction(mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  summarizerFunction(                    outcome_k1,
                                                         outcome_k2,
                                                         outcome_functions,
                                                         cnuk = cnuk, mitj = mitj,
                                                         range_maxi = range_max,
                                                         range_mini = range_min, verbose = verbose,
                                                         memory_manage = memory_manage,
                                                         skinny_summarize_call = skinny_summarize_call,
                                                         SEN=SEN
                  )
                gc()
              }
              gc()
              dataRelatedPair_merge <- readr::read_rds("dataRelatedPair_merge.RDS") # is needed when there are multiple cnu in the loop
            } else { # unoptimized
              temp <- debug_tbl <- try_NA(dataRelatedPair_merge %>% dplyr::filter(.data$cnuRel == cnuk) %>%
                                            sliceFunction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                                            mutateFunction(mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage))

              if (verbose) {
                message("Columns after mutate:")
                print(names(debug_tbl))



                poly_vars_k1 <-     outcome_k1[outcome_functions %in% c("polychorFunction", "ml_polychorFunction")]
                poly_vars_k2 <-     outcome_k2[outcome_functions %in% c("polychorFunction", "ml_polychorFunction")]

                if (length(poly_vars_k1) > 0) {
                  v <- poly_vars_k1[1]
                  v1 <- paste0(poly_vars_k1, "_k1")
                  v2 <- paste0(poly_vars_k2, "_k2")

                  if (all(c(v1, v2, "casecontrol_groupings") %in% names(debug_tbl))) {
                    message(paste(v, "table, cases:"))
                    print(table(
                      debug_tbl[[v1]][debug_tbl$casecontrol_groupings == 1],
                      debug_tbl[[v2]][debug_tbl$casecontrol_groupings == 1],
                      useNA = "ifany"
                    ))

                    message(paste(v, "table, controls:"))
                    print(table(
                      debug_tbl[[v1]][debug_tbl$casecontrol_groupings == 0],
                      debug_tbl[[v2]][debug_tbl$casecontrol_groupings == 0],
                      useNA = "ifany"
                    ))
                  }
                }
              }

              temp <- try_NA(temp   %>%
                               group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                               summarizerFunction(outcome_k1,
                                                  outcome_k2,
                                                  outcome_functions,
                                                  cnuk = cnuk, mitj = mitj,
                                                  range_maxi = range_max,
                                                  range_mini = range_min,
                                                  verbose = verbose, memory_manage = memory_manage,
                                                  skinny_summarize_call = skinny_summarize_call,
                                                  SEN=SEN
                               ))
              # %>%
              # suppressWarnings()
              gc()
            }

          }

          # skip row writing if NA
          if (length(temp) == 1 && is.na(temp)) {
            if (verbose) {
              message(paste0("Skipped writing Temp to disk. Temp was NA"), length(temp))
            }
          } else {
            data.table::fwrite(temp,
                   file = file_path_txt, sep = ",",
                   append = TRUE, row.names = FALSE, col.names = FALSE
            )
            if (memory_manage > 0L) {
              gc()
            }

            # temp get names
            file_names <- names(temp)
            remove(temp)
          }
          # optional loop by grouping variable?
        } # end cnu
        remove(dataRelatedPair_merge)
      } # end mit
      #  if(verbose){
      message(paste(i, addRel_mins[i]))
      #  }
      if (memory_manage > 0L) {
        gc()
      }
    } # end add

    # get the full file to add variable names
    aim1_cors <- utils::read.csv(file_path_txt, header = FALSE)


    names(aim1_cors) <- file_names
    data.table::fwrite(aim1_cors,
           file = file_path_txt, sep = ",",
           append = FALSE, row.names = FALSE, col.names = TRUE
    )
  } # end bin
  if (memory_manage > 0L) {
    gc()
  }
  # clean up after
  remove(dataRelatedPair_merge)
  remove(temp)
  remove(aim1_cors)
  if (cleanup) {
    if (file.exists("dataRelatedPair_merge.RDS")) {
      file.remove("dataRelatedPair_merge.RDS")
    }
    if (file.exists("dataRelatedPair_merge_1.RDS")) {
      file.remove("dataRelatedPair_merge_1.RDS")
    }
    if (file.exists("dataRelatedPair_merge_0.RDS")) {
      file.remove("dataRelatedPair_merge_0.RDS")
    }
  }
  gc()
}
# note: these kin are double entered

