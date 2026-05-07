#------------------------------------------------------------------------------
# not in
`%notin%` <- Negate(`%in%`)



#try return NA

try_NA = function(expr){
  tryCatch(expr,error=function(err) NA)
}


mp_id_to_numeric <- function(x, prefix = "MP") {
  x <- if (!is.character(x)) as.character(x) else x
  y <- sub(paste0("^", prefix), "", x, perl = TRUE)
  suppressWarnings(as.numeric(y))
}
