#------------------------------------------------------------------------------

#' Not-In Operator
#'
#' A convenience operator that returns \code{TRUE} for each element of \code{x}
#' that is \emph{not} present in \code{table}.  It is the negation of \code{\%in\%}.
#'
#' @param x A vector of values to test.
#' @param table A vector of values to test against.
#'
#' @return A logical vector the same length as \code{x}: \code{TRUE} where
#'   \code{x[i]} is not in \code{table}, \code{FALSE} otherwise.
#'
#' @keywords internal
`%notin%` <- Negate(`%in%`)


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
