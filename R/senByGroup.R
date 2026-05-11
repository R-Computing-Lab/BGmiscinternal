#------------------------------------------------------------------------------
# Author Mason Garrison
# Date: 815
# Filename: functions_SENByGroup
# Purpose: this code calcuates the correlation between the outcomes for each kin group, groups by mtdna, cnu, and bins of R

#------------------------------------------------------------------------------

## First are the helper functions
options(scipen = 10, digits = 11)
## passes subset if slice_1000 is true, otherwise passes entire thing

source("functions/14CorrelateOutcomesByGroup_helperfunctions.R")

# here's the mega function
SENByGroup <- function(df_foldername = "longevity_skinny_matpat",
                       binwidth_num = c(.1, .05),
                       binwidth_cha = c("10", "05"),
                       kin_degree_max = 12,
                       kin_degree_min = 0,
                       max_kin_per_bin = 8.7 * 10^9, # 10^7 is 10 million
                       cnu = c(1, 0),
                       mit = c(1, 0),
                       doubleentered = TRUE,
                       slice_1000 = TRUE,
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
                         "USA_flag_10",
                         "USA_flag_10"
                       ),
                       outcome_functions = c(
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
                       age_filter = FALSE, min_age = 0) {
  # potential options
  # expand.grid on unique values and filter via a loop
  ## make sure to use verbose to give folks a sense of what they;re asking for
  ## include a check to estimate time?
  # groupby
  options(scipen = 10, digits = 11)

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
          dataRelatedPair_merge <- NA
          next
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
                dplyr::filter(age_k1 >= min_age | age_k2 >= min_age)
            } else if ("age" %in% names(dataRelatedPair_merge)) {
              dataRelatedPair_merge <- dataRelatedPair_merge %>%
                drop_na(age) %>%
                dplyr::filter(age >= min_age)
            }
            if (verbose == TRUE) {
              message(paste(
                i, addRel_mins[i],
                "Filtering people who didn't live to",
                min_age, " ", current_rows - nrow(dataRelatedPair_merge), "removed"
              ))
            }
          }



          # possible optimization
          dataRelatedPair_merge <- dataRelatedPair_merge %>%
            dplyr::mutate(dplyr::across(where(is.numeric), round, 3))
        }

        if (length(dataRelatedPair_merge) == 0 || nrow(dataRelatedPair_merge) > max_kin_per_bin || (length(dataRelatedPair_merge) == 1 && is.na(dataRelatedPair_merge))) {
          # this loop skips the file if it is too big, otherwise the loop fails
          # and gives unhappy warning and fails without cleaning up
          message(paste(
            i, addRel_mins[i],
            "was skipped because it was ",
            nrow(dataRelatedPair_merge),
            "which is bigger than ", max_kin_per_bin
          ))
          remove(dataRelatedPair_merge)
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
            } else { # if IDs not there, this is a painlessway to reduce the size of the matrix
              dxlist <- c( # "ID1","ID2", #intentional ordering
                "addRel", # "mitRel",
                "cnuRel",
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k2")],
                names(dataRelatedPair_merge)[endsWith(names(dataRelatedPair_merge), "_k1")]
              )
            }
            if (verbose) {
              print(dxlist)
            }
            # remove addRel
            if ("addRel" %notin% names(dataRelatedPair_merge)) {
              dxlist <- dxlist[!dxlist %in% "addRel"]
            }

            if ("cnuRel" %notin% names(dataRelatedPair_merge)) {
              dxlist <- dxlist[!dxlist %in% "cnuRel"]
            }


            dataRelatedPair_merge <- rbindlist(
              list(
                dataRelatedPair_merge,
                dataRelatedPair_merge[, ..dxlist]
              ),
              use.names = FALSE
            )
            gc()
            remove(dxlist)
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
            next
          } else {
            if (verbose) {
              #   message(head(dataRelatedPair_merge))
            }
            #  message("here")
            # sometimes the nested functions can't find the looping index
            cnuk <- cnu[k]
            mitj <- mit[j]

            temp <- try_NA(dataRelatedPair_merge %>% dplyr::filter(cnuRel == cnuk) %>%
                             sliceFuction(slice_1000 = slice_1000, memory_manage = memory_manage) %>%
                             mutateFunction(mutate_vars = mutate_vars, verbose = verbose, memory_manage = memory_manage) %>%
                             #     group_by(ID1) %>% dplyr::mutate(
                             #     unique_ID2s = dplyr::n_distinct(ID2)
                             #      ) %>% ungroup() %>% group_by(ID2) %>% dplyr::mutate(
                             #        unique_ID1s = dplyr::n_distinct(ID1)
                             #      ) %>% ungroup()  %>%
                             group_byFunction(grouping_vars, verbose = verbose, memory_manage = memory_manage) %>%
                             summarizerFunction(outcome_vars, outcome_functions,
                                                cnuk = cnuk, mitj = mitj,
                                                range_maxi = range_max,
                                                range_mini = range_min,
                                                verbose = verbose, memory_manage = memory_manage,
                                                skinny_summarize_call = skinny_summarize_call,
                                                SEN=TRUE
                             ))
            # %>%
            # suppressWarnings()
            gc()
          }

          # }

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
    } # end add

    # get the full file to add variable names
    aim1_cors <- read.csv(file_path_txt, header = FALSE)


    names(aim1_cors) <- file_names
    fwrite(aim1_cors,
           file = file_path_txt, sep = ",",
           append = FALSE, row.names = FALSE, col.names = TRUE
    )
  } # end bin

  # clean up after
  remove(dataRelatedPair_merge)
  remove(temp)
  remove(aim1_cors)

}
# note: these kin are double entered
estimate_icc_latent_from_dyadic <- function(tbl, outcome_var, method = c("latent", "mean", "glmer","lmer")) {
  method <- match.arg(method)

  outcome_k1 <- paste0(outcome_var, "_k1")
  outcome_k2 <- paste0(outcome_var, "_k2")

  df_long <- tibble::tibble(
    ID = c(tbl$ID1, tbl$ID2),
    outcome = c(tbl[[outcome_k1]], tbl[[outcome_k2]])
  ) %>%
    dplyr::filter(!is.na(outcome))


}
estimate_icc_latent_from_dyadic <- function(tbl, outcome_var,
                                            method = c("latent", "mean", "lmer","glmer"),
                                            binary = TRUE) {
  method <- match.arg(method)

  outcome_k1 <- paste0(outcome_var, "_k1")
  outcome_k2 <- paste0(outcome_var, "_k2")

  df_long <- tibble::tibble(
    ID = c(tbl$ID1, tbl$ID2),
    outcome = c(tbl[[outcome_k1]], tbl[[outcome_k2]])
  ) %>%
    dplyr::filter(!is.na(outcome))

  if (method == "mean") {
    # ICC estimated as ratio of between to total variance of per-ID means
    df_id <- df_long %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(mu = mean(outcome), .groups = "drop")

    between_var <- var(df_id$mu)
    total_var <- var(df_long$outcome)

    icc <- between_var / total_var
    return(icc)

  } else if (method == "latent") {
    # For binary: ICC on logistic latent scale; for continuous: standard ICC
    df_long <- df_long %>%
      dplyr::group_by(ID) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup()

    person_means <- df_long %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(mean_outcome = mean(outcome), .groups = "drop")

    between_var <- var(person_means$mean_outcome)

    if (binary) {
      p_bar <- mean(df_long$outcome)
      total_var <- p_bar * (1 - p_bar)
    } else {
      total_var <- var(df_long$outcome)
    }

    icc <- between_var / total_var
    return(icc)

  } else if (method=="lmer"){

    icc_model <- lme4::lmer(outcome ~ 1 + (1 | ID), data = df_long)

    # Extract variance components
    var_components <- as.data.frame(VarCorr(icc_model))
    sigma_b <- var_components$vcov[1]
    sigma_e <- attr(VarCorr(icc_model), "sc")^2

    icc <- sigma_b / (sigma_b + sigma_e)
    return(icc)

  } else if (method=="psych"){
    library(psych)
    icc_data <- df_long %>%
      group_by(ID) %>%
      dplyr::mutate(obs_num = row_number()) %>%
      ungroup() %>%
      tidyr::pivot_wider(names_from = obs_num, values_from = outcome)

    # compute ICC (only works if enough repeated measures per person)
    icc <- psych::ICC(icc_data[,-1])  # drop ID column
    return(icc)
  } else if (method == "glmer") {
    if (binary==FALSE) {
      stop("The 'glmer' method is only valid for binary outcomes (binary = TRUE).")
    }

    df_long <- df_long %>%
      dplyr::group_by(ID) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup()

    model <- lme4::glmer(outcome ~ 1 + (1 | ID), data = df_long, family = binomial)
    vc <- as.data.frame(VarCorr(model))$vcov[1]
    icc <- vc / (vc + (pi^2 / 3))  # Latent scale variance for logistic link
    return(icc)
  }
}
