#------------------------------------------------------------------------------


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
  tryCatch(expr, error = function(err){return(NA)})
}
#' Not-In Operator
#'
#' A convenience operator that returns \code{TRUE} for each element of \code{x}
#' that is \emph{not} present in \code{table}.  It is the negation of \code{\%in\%}.
#'
#'
#' @return A logical vector the same length as \code{x}: \code{TRUE} where
#'   \code{x[i]} is not in \code{table}, \code{FALSE} otherwise.
#'
#' @keywords internal
`%notin%` <- Negate(`%in%`)



#' Convert MP ID to Numeric
#'
#' Strips a leading prefix (default \code{"MP"}) from character IDs and coerces
#' the remainder to \code{numeric}.
#'
#' @param x A character vector (or an object coercible to character) of IDs.
#' @param prefix A single string giving the prefix to remove.  Default is
#'   \code{"MP"}.
#'
#' @return A numeric vector the same length as \code{x}.  Elements that cannot
#'   be coerced after prefix removal are returned as \code{NA_real_}.
#'
#' @export
#' @examples
#' mp_id_to_numeric(c("MP001", "MP042", "999"))
#' mp_id_to_numeric(c("ID10", "ID20"), prefix = "ID")
mp_id_to_numeric <- function(x, prefix = "MP") {
  x <- if (!is.character(x)) as.character(x) else x
  y <- sub(paste0("^", prefix), "", x, perl = TRUE)
  suppressWarnings(as.numeric(y))
}


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


