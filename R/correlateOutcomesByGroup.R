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
                                     outcome_vars = c(
                                       "male", "age",
                                       "USA_flag_10",
                                       "USA_flag_15",
                                       "USA_flag_10",
                                       "USA_flag_10",
                                       "USA_flag_15",
                                       "USA_flag_15"
                                     ),
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
                drop_na(age_k1, age_k2) %>%
                filter(age_k1 >= min_age & age_k2 >= min_age)
            } else if ("age" %in% names(dataRelatedPair_merge)) {
              dataRelatedPair_merge <- dataRelatedPair_merge %>%
                drop_na(age) %>%
                filter(age >= min_age)
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
              dplyr::mutate(across(where(is.numeric), round, 3))

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
              print(names(dataRelatedPair_merge[, ..dxlist_main]))
              print(names(dataRelatedPair_merge[, ..dxlist]))
            }
            dataRelatedPair_merge <- rbindlist(
              list(
                dataRelatedPair_merge[, ..dxlist_main],
                dataRelatedPair_merge[, ..dxlist]
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
          current_nrows <- nrow(dataRelatedPair_merge %>% dplyr::filter(cnuRel == cnu[k]))

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
              require(tidyft)
              # this writes each iteration to disk the idea being that you loop thru all the groups
              ## not great solution but... it does let you use mutate once

              dataRelatedPair_merge %>%
                mutateFunction(
                  mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage
                ) %>%
                convert_to_integer(memory_manage = memory_manage) %>%
                write_rds("dataRelatedPair_merge.RDS")

              dataRelatedPair_merge <- tidyft::setDT(dataRelatedPair_merge) %>% tidyft::filter(cnuRel == cnuk)
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
                  dataRelatedPair_merge <- read_rds("dataRelatedPair_merge.RDS") %>%
                    tidyft::setDT() %>%
                    tidyft::filter(cnuRel == cnuk)
                  gc()
                } else {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- read_rds("dataRelatedPair_merge.RDS") %>%
                    tidyft::setDT(dataRelatedPair_merge) %>%
                    tidyft::filter(.data[[grouping_vars]] == grouping_loop[g], cnuRel == cnuk)
                  gc()
                }

                temp <- dataRelatedPair_merge %>%
                  sliceFuction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                  #    mutateFunction(mutate_vars=mutate_vars,verbose=verbose, memory_manage = memory_manage)%>%
                  group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  summarizerFunction(outcome_vars,
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
              dataRelatedPair_merge <- read_rds("dataRelatedPair_merge.RDS")
            } else if (memory_manage > 0L) { # clean up the this to be one variable with 3 levels, character is fine, 0L
              gc()

              # this writes each iteration to disk the idea being that you loop thru all the groups
              ## not great solution but... it does let you use mutate once
              # print("353")
              dataRelatedPair_merge %>%
                # mutateFunction(
                #    mutate_vars = mutate_vars, verbose = verbose,  memory_manage= memory_manage
                # ) %>% convert_to_integer(memory_manage =  memory_manage) %>%
                write_rds("dataRelatedPair_merge.RDS")

              dataRelatedPair_merge %>%
                mutateFunction(
                  mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage
                ) %>%
                dplyr::filter(cnuRel == cnuk) %>%
                convert_to_integer(memory_manage = memory_manage) %>%
                select(-cnuRel) %>%
                write_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS"))


              dataRelatedPair_merge <- dataRelatedPair_merge %>% dplyr::filter(cnuRel == cnuk)
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
                  dataRelatedPair_merge <- read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) # %>% # this can be managed
                  #  dplyr::filter(cnuRel == cnuk)
                  gc()
                } else if (FALSE) { # &exists("grouping_loop_grid")){
                  # if there's a grid
                  ## not complete
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) %>%
                    dplyr::filter(.data[[grouping_vars[1]]] == grouping_loop_grid[1, g], .data[[grouping_vars[2]]] == grouping_loop_grid[2, g])
                  gc()
                } else {
                  remove(dataRelatedPair_merge)
                  gc()
                  dataRelatedPair_merge <- read_rds(paste0("dataRelatedPair_merge_", cnuk, ".RDS")) %>%
                    dplyr::filter(.data[[grouping_vars]] == grouping_loop[g])
                  gc()
                }
                temp <- dataRelatedPair_merge %>%
                  sliceFuction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                  mutateFunction(mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                  summarizerFunction(outcome_vars,
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
              dataRelatedPair_merge <- read_rds("dataRelatedPair_merge.RDS") # is needed when there are multiple cnu in the loop
            } else { # unoptimized
              temp <- debug_tbl <- try_NA(dataRelatedPair_merge %>% dplyr::filter(cnuRel == cnuk) %>%
                                            sliceFuction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                                            mutateFunction(mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage))

              if (verbose) {
                message("Columns after mutate:")
                print(names(debug_tbl))



                poly_vars <- outcome_vars[outcome_functions %in% c("polychorFunction", "ml_polychorFunction")]

                if (length(poly_vars) > 0) {
                  v <- poly_vars[1]
                  v1 <- paste0(v, "_k1")
                  v2 <- paste0(v, "_k2")

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
                               summarizerFunction(outcome_vars, outcome_functions,
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
            fwrite(temp,
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
    aim1_cors <- read.csv(file_path_txt, header = FALSE)


    names(aim1_cors) <- file_names
    fwrite(aim1_cors,
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

