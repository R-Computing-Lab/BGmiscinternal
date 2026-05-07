#------------------------------------------------------------------------------
# Author Mason Garrison
# Date: 815
# Filename: functions_CorrelateOutcomesByGroup
# Purpose: this code calcuates the correlation between the outcomes for each kin group, groups by mtdna, cnu, and bins of R

#------------------------------------------------------------------------------

## First are the helper functions
options(scipen = 10, digits = 11)


# not in
`%notin%` <- Negate(`%in%`)

# convert numeric to integers

convert_to_integer <- function(tbl,memory_manage = 0L){
  if(memory_manage<1){
    tbl %>%
      mutate(across(where(is.numeric), ~ if_else(. == floor(.), as.integer(.),.)))
    #   require(varhandle)
  }else{
    require(data.table)
    # convert to datatable
    tbl <- data.table::as.data.table(tbl)
    # convert numeric to integer where possible
    for(col in names(tbl)) {
      if(is.numeric(tbl[[col]])){
        # if intger the same for both, then convert in place
        if(mean(floor(unique(tbl[[col]]))==unique(tbl[[col]]),na.rm=TRUE)){
          tbl[,(col) := as.integer(tbl[[col]])]
        }
      }
    }
    as_tibble(tbl)
  }
}

#try return NA

try_NA = function(expr){
  tryCatch(expr,error=function(err) NA)
}

## passes subset if slice_1000 is true, otherwise passes entire thing
sliceFuction <- function(tbl, slice_1000, memory_manage = 0L) {
  if (slice_1000==TRUE) {
    tbl %>% slice(1:1000)
  } else {
    tbl
  }
}

# passes thru multiple grouping variables. right now I'm using groupby, but this could be replaced by expand grid.
# it does allow for multiple groups
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
mutateFunction <- function(tbl, mutate_vars, verbose = FALSE,  memory_manage = 0L ) {
  # subfunctions

  ## case control

  mutateCaseControl <- function(tbl, memory_manage=0L){
    if ("casecontrol_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
      }
      tbl %>% # match
        dplyr::mutate(
          casecontrol_groupings = case_when(
            CaseControl_k1 == 1 & Has_Dementia_k1 == 1 ~ 1L,
            CaseControl_k1 == 1 & Has_Dementia_k1 == 0 ~ 0L,
            TRUE ~ NA_integer_
          )
        )
    }
  }


  ## gender
  mutateGender <- function(tbl, memory_manage=0L){
    if ("gender_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            gender_groupings = case_when(
              male_k1 == 1 & male_k2 == 1 ~ 2L,
              male_k1 == 0 & male_k2 == 0 ~ 0L,
              male_k1 == 1 & male_k2 == 0 ~ 1L,
              male_k1 == 0 & male_k2 == 1 ~ 1L,
              TRUE ~ NA_integer_
            )
          )   %>% select(-c(male_k2))
      } else{
        tbl %>% # match
          dplyr::mutate(
            gender_groupings = case_when(
              male_k1 == 1 & male_k2 == 1 ~ 2L,
              male_k1 == 0 & male_k2 == 0 ~ 0L,
              male_k1 == 1 & male_k2 == 0 ~ 1L,
              male_k1 == 0 & male_k2 == 1 ~ 1L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  ## linkage type
  mutateLinkage <- function(tbl, memory_manage=0L){
    if ("linkagetype" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )) %>% select(-c(patID_k1,patID_k2,matID_k1,matID_k2)) %>%
          dplyr::mutate(  linkagetype = case_when(
            same_matID == 1 & same_patID == 1 ~ 11L,
            same_matID == 0 & same_patID == 0 ~ 00L,
            same_matID == 1 & same_patID == 0 ~ 10L,
            same_matID == 0 & same_patID == 1 ~ 01L,
            TRUE ~ NA_integer_
          )) %>% select(-c(same_matID, same_patID))
      } else{
        tbl %>% # match
          dplyr::mutate(
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            linkagetype = case_when(
              same_matID == 1 & same_patID == 1 ~ 11L,
              same_matID == 0 & same_patID == 0 ~ 00L,
              same_matID == 1 & same_patID == 0 ~ 10L,
              same_matID == 0 & same_patID == 1 ~ 01L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }

  mutateLinkage_same_patID <- function(tbl, memory_manage=0L){
    if ("same_patID" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            #  same_matID = case_when(
            #    matID_k1 == matID_k2 ~ 1L,
            #    matID_k1 != matID_k2 ~ 0L,
            #    TRUE ~ NA_integer_
            #  ),
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )) %>% select(-c(patID_k1,patID_k2,matID_k1,matID_k2))
      } else{
        tbl %>% # match
          dplyr::mutate(
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ))
      }
    }
  }
  mutateLinkage_same_matID <- function(tbl, memory_manage=0L){
    if ("same_matID" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            #  same_matID = case_when(
            #    matID_k1 == matID_k2 ~ 1L,
            #    matID_k1 != matID_k2 ~ 0L,
            #    TRUE ~ NA_integer_
            #  ),
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )) %>% select(-c(patID_k1,patID_k2,matID_k1,matID_k2))
      } else{
        tbl %>% # match
          dplyr::mutate(
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ))
      }
    }
  }

  mutateLinkage_any <- function(tbl, memory_manage=0L){
    if ("linkagetype" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            )) %>% select(-c(patID_k1,patID_k2,matID_k1,matID_k2)) %>%
          dplyr::mutate(  linkagetype = case_when(
            same_matID == 1 & same_patID == 1 ~ 11L,
            same_matID == 0 & same_patID == 0 ~ 00L,
            same_matID == 1 & same_patID == 0 ~ 11L,
            same_matID == 0 & same_patID == 1 ~ 11L,
            TRUE ~ NA_integer_
          )) %>% select(-c(same_matID, same_patID))
      } else{
        tbl %>% # match
          dplyr::mutate(
            same_matID = case_when(
              matID_k1 == matID_k2 ~ 1L,
              matID_k1 != matID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            same_patID = case_when(
              patID_k1 == patID_k2 ~ 1L,
              patID_k1 != patID_k2 ~ 0L,
              TRUE ~ NA_integer_
            ),
            linkagetype = case_when(
              same_matID == 1 & same_patID == 1 ~ 11L,
              same_matID == 0 & same_patID == 0 ~ 00L,
              same_matID == 1 & same_patID == 0 ~ 11L,
              same_matID == 0 & same_patID == 1 ~ 11L,
              TRUE ~ NA_integer_
            )
          )
      }
    }
  }
  mutateCohort_19 <- function(tbl, memory_manage=0L){
    if ("cohort_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            )) %>% select(-c(BYr_k1,BYr_k2)) %>%
          dplyr::mutate(
            cohort_groupings = case_when(cohort_k1 < cohort_k2 ~ as.integer(cohort_k1*100+cohort_k2),
                                         cohort_k1 >= cohort_k2 ~ as.integer(cohort_k2*100+cohort_k1),
                                         TRUE ~ NA_integer_)
          ) %>% select(-c(cohort_k1,cohort_k2))
      }else{ tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              cohort_k1 < cohort_k2, paste(cohort_k1, cohort_k2, sep = "_"),
              paste(cohort_k2, cohort_k1, sep = "_")
            )
          )
      }
    }
  }
  mutateCohort_19flat <- function(tbl, memory_manage=0L){
    if ("cohort_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            )) %>% select(-c(BYr_k1,BYr_k2)) %>%
          dplyr::mutate(
            cohort_groupings = case_when(cohort_k1 == 19 & cohort_k2 == 19 ~ as.integer(1919),
                                         cohort_k1 < 19 | cohort_k2 <  19 ~ as.integer(1818),
                                         TRUE ~ NA_integer_)
          ) %>% select(-c(cohort_k1,cohort_k2))
      }else{ tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 ~ 19L,
              TRUE ~ NA_integer_
            ), cohort_groupings = case_when(cohort_k1 == 19 & cohort_k2 == 19 ~ as.integer(1919),
                                            cohort_k1 < 19 | cohort_k2 <  19 ~ as.integer(1818),
                                            TRUE ~ NA_integer_)
          )
      }
    }
  }
  # cohort_century
  mutateCohortMatch <- function(tbl, memory_manage=0L){
    if ("cohort_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 >= 1600 & BYr_k1 < 1700 ~ 16L,
              BYr_k1 >= 1700 & BYr_k1 < 1800 ~ 17L,
              BYr_k1 >= 1800 & BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 & BYr_k1 < 2000 ~ 19L,
              BYr_k1 >= 2000 & BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 >= 1600 & BYr_k2 < 1700 ~ 16L,
              BYr_k2 >= 1700 & BYr_k2 < 1800 ~ 17L,
              BYr_k2 >= 1800 & BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 & BYr_k2 < 2000 ~ 19L,
              BYr_k2 >= 2000 & BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            )) %>% select(-c(BYr_k1,BYr_k2)) %>%
          dplyr::mutate(
            cohort_groupings = case_when(cohort_k1 == cohort_k2 ~ as.integer(1),
                                         cohort_k1 != cohort_k2 ~ as.integer(0),
                                         TRUE ~ NA_integer_)
          ) %>% select(-c(cohort_k1,cohort_k2))
      }else{ tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 >= 1600 & BYr_k1 < 1700 ~ 16L,
              BYr_k1 >= 1700 & BYr_k1 < 1800 ~ 17L,
              BYr_k1 >= 1800 & BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 & BYr_k1 < 2000 ~ 19L,
              BYr_k1 >= 2000 & BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 >= 1600 & BYr_k2 < 1700 ~ 16L,
              BYr_k2 >= 1700 & BYr_k2 < 1800 ~ 17L,
              BYr_k2 >= 1800 & BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 & BYr_k2 < 2000 ~ 19L,
              BYr_k2 >= 2000 & BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              cohort_k1 == cohort_k2, as.integer(1),
              as.integer(0)
            )
          )
      }
    }
  }

  # cohort_century
  mutateCohort <- function(tbl, memory_manage=0L){
    if ("cohort_groupings" %in% names(tbl)){
      tbl # skip in outcome vare already present
    } else{
      if (memory_manage > 0L) {
        gc()
        tbl %>% # match
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 >= 1600 & BYr_k1 < 1700 ~ 16L,
              BYr_k1 >= 1700 & BYr_k1 < 1800 ~ 17L,
              BYr_k1 >= 1800 & BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 & BYr_k1 < 2000 ~ 19L,
              BYr_k1 >= 2000 & BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 >= 1600 & BYr_k2 < 1700 ~ 16L,
              BYr_k2 >= 1700 & BYr_k2 < 1800 ~ 17L,
              BYr_k2 >= 1800 & BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 & BYr_k2 < 2000 ~ 19L,
              BYr_k2 >= 2000 & BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            )) %>% select(-c(BYr_k1,BYr_k2)) %>%
          dplyr::mutate(
            cohort_groupings = case_when(cohort_k1 < cohort_k2 ~ as.integer(cohort_k1*100+cohort_k2),
                                         cohort_k1 >= cohort_k2 ~ as.integer(cohort_k2*100+cohort_k1),
                                         TRUE ~ NA_integer_)
          ) %>% select(-c(cohort_k1,cohort_k2))
      }else{ tbl %>% # match
          # birthcohort
          dplyr::mutate(
            cohort_k1 = case_when(
              BYr_k1 >= 1600 & BYr_k1 < 1700 ~ 16L,
              BYr_k1 >= 1700 & BYr_k1 < 1800 ~ 17L,
              BYr_k1 >= 1800 & BYr_k1 < 1900 ~ 18L,
              BYr_k1 >= 1900 & BYr_k1 < 2000 ~ 19L,
              BYr_k1 >= 2000 & BYr_k1 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ),
            cohort_k2 = case_when(
              BYr_k2 >= 1600 & BYr_k2 < 1700 ~ 16L,
              BYr_k2 >= 1700 & BYr_k2 < 1800 ~ 17L,
              BYr_k2 >= 1800 & BYr_k2 < 1900 ~ 18L,
              BYr_k2 >= 1900 & BYr_k2 < 2000 ~ 19L,
              BYr_k2 >= 2000 & BYr_k2 < 2100 ~ 20L,
              TRUE ~ NA_integer_
            ), cohort_groupings = ifelse(
              cohort_k1 < cohort_k2, paste(cohort_k1, cohort_k2, sep = "_"),
              paste(cohort_k2, cohort_k1, sep = "_")
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
  } else if (mutate_vars == "gender_groupings_linkagetype"|mutate_vars == "gender_linkagetype") {
    tbl %>% mutateLinkage(memory_manage=memory_manage) %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mutate_vars == "linkagetype") {
    tbl %>% # match
      mutateLinkage(memory_manage=memory_manage)
  } else if (mutate_vars == "linkage_any"){
    tbl %>% # match
      mutateLinkage_any(memory_manage=memory_manage)
  } else if (mutate_vars == "same_patID") {
    tbl %>% # match
      mutateLinkage_same_patID(memory_manage=memory_manage)
  } else if (mutate_vars == "same_matID") {
    tbl %>% # match
      mutateLinkage_same_matID(memory_manage=memory_manage)
  } else if (mutate_vars == "gender_groupings") {
    tbl %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mutate_vars == "casecontrol_groupings") {
    tbl %>%
      mutateCaseControl(memory_manage = memory_manage)
  } else if (mutate_vars == "cohort_groupings_19") {
    tbl %>% mutateCohort_19(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_groupings_19flat") {
    tbl %>% mutateCohort_19flat(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_groupings_match") {
    tbl %>% mutateCohortMatch(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_groupings") {
    tbl %>% mutateCohort(memory_manage=memory_manage)
  } else if (mutate_vars == "gender_cohort_groupings_19") {
    tbl %>% mutateCohort_19(memory_manage=memory_manage) %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mutate_vars == "gender_cohort_groupings") {
    tbl %>% mutateCohort(memory_manage=memory_manage) %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_linkagetype_19") {
    tbl %>% mutateLinkage(memory_manage=memory_manage) %>%
      mutateCohort_19(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_linkagetype") {
    tbl %>% mutateLinkage(memory_manage=memory_manage) %>%
      mutateCohort(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_gender_linkagetype_19") {
    tbl %>% mutateLinkage(memory_manage=memory_manage) %>%
      mutateCohort_19(memory_manage=memory_manage) %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mutate_vars == "cohort_gender_linkagetype") {
    tbl %>% mutateLinkage(memory_manage=memory_manage) %>%
      mutateCohort(memory_manage=memory_manage) %>%
      mutateGender(memory_manage=memory_manage)
  } else if (mode(mutate_vars) == "character") {
    mutate_call <- paste0("tbl %>% mutate(", mutate_vars, ")")
    eval(parse(text = mutate_call))
  } else {
    tbl
  }
}





# creates the group bys
## I tried A TON OF THINGS, but the most predictably behaving version...
## was to create a large string via mapply
summarizerFunction <- function(tbl,
                               outcome_k1,
                               outcome_k2, outcome_functions, mitj, cnuk,
                               range_maxi = range_max,
                               range_mini, verbose = FALSE, memory_manage = 0L,skinny_summarize_call=TRUE,
                               SEN=FALSE
) {
  if (length(outcome_k1) != length(outcome_functions)) {
    stop("The vectors of function names and variables must be the same length")
  }
  if (length(outcome_k1) != length(outcome_k2)) {
    stop("The vectors of outcome names must be the same length")
  }

  if(SEN==FALSE){
    if(skinny_summarize_call==TRUE){
      if(memory_manage>1|| !("addRel" %in% names(tbl))) {
        summarize_call <- "tbl %>% summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      } else{
        summarize_call <- "tbl %>% summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
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
      if(memory_manage>1|| !("addRel" %in% names(tbl))) {
        summarize_call <- "tbl %>% summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = n_distinct(c(ID2,ID1)),
	    addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"
      } else{
        summarize_call <- "tbl %>% summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = n_distinct(c(ID2,ID1)),
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
   summarize( n_pairs = n()/(1+doubleentered), # if double entered this value is 2, and if not double entered this value is 1)
      unique_n = n_distinct(c(ID1, ID2)),
      avg_rows_per_id = n() / n_distinct(c(ID1, ID2)),
      avg_dyads_per_id = (n() / (1 + doubleentered)) / n_distinct(c(ID1, ID2)),
    #  icc_by_id = estimate_icc_latent_from_dyadic(outcome_var='USA_flag_10', method = 'mean'),
	    addRel_min = try_NA(range_mini),
      addRel_max = try_NA(range_maxi),
      mtdna = try_NA(mitj),
      cnu = try_NA(cnuk),"

  }


  summarize_parts <- mapply(function(var_k1,
                                     var_k2,
                                     fun) {
    if (fun == "polychorFunction") {
      paste0(var_k1, "_", fun, " = list(try_NA(polychor(",
             var_k1, "_k1,",
             var_k2, "_k2, std.err=TRUE)) %>%
                 {list(rho = try_NA(.$rho),
                       se = sqrt(try_NA(.$var)),
                       chisq = try_NA(.$chisq),
                       df = try_NA(.$df))})")
    } else if (fun == "ml_polychorFunction") {
      paste0(var_k1, "_", fun, " = try_NA(polychor(", var_k1, "_k1,", var_k2, "_k2,ML=TRUE))")
    } else if (fun == "relriskFunction") {
      paste0(var_k1, "_", fun, " = list(try_NA(relriskFunction(",
             var_k1, "_k1, ",
             var_k2, "_k2)) %>%
             {list(rr = try_NA(.[1]),
                   LL = try_NA(.[2]),
                   UL = try_NA(.[3]))})")
    } else if (fun %in% c("phi_both", "rr_exposed_cases", "rr_a")){
      paste0(var_k1, "_", fun, "= try_NA(sum(", var_k1, "_k1 == 1 & ",var_k2, "_k2 == 1, na.rm=TRUE))")
    } else if (fun %in% c("phi_none", "rr_unexposed_noncases", "rr_d")) {
      paste0(var_k1, "_", fun, "= try_NA(sum(", var_k1, "_k1 == 0 & ",var_k2, "_k2 == 0, na.rm=TRUE))")
    } else if (fun %in% c("phi_one", "rr_discordant", "rr_b_plus_c")) {
      paste0(var_k1, "_", fun, "= try_NA(sum((", var_k1, "_k1 == 1 & ",var_k2, "_k2 == 0) | (", var_k1, "_k1 == 0 & ",var_k2, "_k2 == 1), na.rm = TRUE))")
    } else if (fun %in% c("phi_k1_yes_k2_no", "rr_exposed_noncases", "rr_b")) {
      paste0(var_k1, "_", fun, " = try_NA(sum(", var_k1, "_k1 == 1 & ", var_k2, "_k2 == 0, na.rm = TRUE))")

    } else if (fun %in% c("phi_k1_no_k2_yes", "rr_unexposed_cases", "rr_c")) {
      paste0(var_k1, "_", fun, " = try_NA(sum(", var_k1, "_k1 == 0 & ", var_k2, "_k2 == 1, na.rm = TRUE))")
    } else if (fun == "phi_est"){
      paste0(var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
             var_k1, "_phi_both, f01= ",
             var_k1, "_phi_one, f10=",
             var_k1, "_phi_one, f00 = ",
             var_k1, "_phi_none))[1]")
    } else if (fun == "phi_se"){
      paste0(var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
             var_k1, "_phi_both, f01= ",
             var_k1, "_phi_one,f10=",
             var_k1, "_phi_one,f00 = ",
             var_k1, "_phi_none))[2]")
    } else if (fun == "phi_LL"){
      paste0(var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
             var_k1, "_phi_both, f01= ",
             var_k1, "_phi_one,f10=",
             var_k1, "_phi_one,f00 = ",
             var_k1, "_phi_none))[3]")
    } else if (fun == "phi_UL"){
      paste0(var_k1, "_", fun, " = try_NA(ci.phi(.05,f11= ",
             var_k1, "_phi_both, f01= ",
             var_k1, "_phi_one,f10=",
             var_k1, "_phi_one,f00 = ",
             var_k1, "_phi_none))[4]")
    } else if (fun == "phi_ci"){
      paste0(var_k1, "_", fun, " = list(try_NA(ci.phi(.05,f11= ",
             var_k1, "_phi_both, f01= ",
             var_k1, "_phi_one,f10=",
             var_k1, "_phi_one,f00 = ",
             var_k1, "_phi_none))  %>%
             {list(phi = try_NA(.$Estimate),
                   se = try_NA(.$SE),
                   LL = try_NA(.$LL),
                   UL = try_NA(.$UL))})")
    } else if (fun == "rr_risk_exposed") {
      paste0(var_k1, "_", fun_k1, " = try_NA(ifelse((",
             var_k1, "_rr_a + ", var_k1, "_rr_b) > 0, ",
             var_k1, "_rr_a / (", var_k1, "_rr_a + ", var_k1, "_rr_b), ",
             "NA_real_))")

    } else if (fun == "rr_risk_unexposed") {
      paste0(var_k1, "_", fun, " = try_NA(ifelse((",
             var_k1, "_rr_c + ", var_k1, "_rr_d) > 0, ",
             var_k1, "_rr_c / (", var_k1, "_rr_c + ", var_k1, "_rr_d), ",
             "NA_real_))")

    } else if (fun == "rr_est") {
      paste0(var_k1, "_", fun, " = try_NA(ifelse(",
             var_k1, "_rr_risk_unexposed > 0, ",
             var_k1, "_rr_risk_exposed / ", var, "_rr_risk_unexposed, ",
             "NA_real_))")
    }else{
      paste0(var_k1, "_", fun, " = try_NA(", fun, "(", var_k1, "))")
    }
  }, outcome_k1, outcome_k2,

  outcome_functions, SIMPLIFY = FALSE)
  # if unnesting variable is needed
  unnest_parts <- c()

  if ("polychorFunction" %in% outcome_functions) {
    # get vars to unlist
    # doesn't like when there are multiple polychors
    var_polychor_unnest <- outcome_vars[outcome_functions == "polychorFunction"]
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


  summarize_call <- paste0(summarize_call, paste(summarize_parts, collapse = ", "), ")", unnest_call)
  if (verbose) {
    message(summarize_call)
  }

  # evaluate the constructed function call and return result
  eval(parse(text = summarize_call))
}


meanFunction <- function(x) {
  mean(x, na.rm = TRUE)
}


sdFunction <- function(x) {
  sd(x, na.rm = TRUE)
}

q25Function <- function(x) {
  quantile(x, na.rm = TRUE,probs = .25)
}
q75Function <- function(x) {
  quantile(x, na.rm = TRUE,probs = .75)
}
q50Function <- function(x) {
  quantile(x, na.rm = TRUE,probs = .5)
}

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


make_input_file <- function(data_path,df_foldername,binwidth_cha,mit,range_min,range_max){
  paste0(data_path,"data/", df_foldername, "_", binwidth_cha, "/df_mt", mit, "_r", range_min, "-r", range_max, ".csv")

}


read_kinbin <- function(input_file, drop_variables = c("mitRel"), verbose = FALSE) {
  dataRelatedPair_merge <- try_NA(fread(input_file,
                                        header = TRUE,
                                        #  drop vars to slim
                                        drop = drop_variables
  )) %>%
    #   mutate(addRel = round(addRel, digits = 4)) %>% # can be dropped?
    suppressWarnings()

  if (verbose == TRUE) {
    message(paste0(input_file, "had ", nrow(dataRelatedPair_merge), " rows"))
  }

  return(dataRelatedPair_merge)
}
