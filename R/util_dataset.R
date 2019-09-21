# Log an operation
logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)

  obj
}

# Conditional mutate
# mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
#   condition <- eval(substitute(condition), .data, envir)
#   .data[condition, ] <- .data[condition, ] %>% mutate(...)
#   .data
# }

# Check if the object has certain fields
has_mobster_fits = function(x)
{
  !all(is.null(x$fit_MOBSTER))
}

has_viber_fits = function(x)
{
  !all(is.null(x$fit_VIBER))
}

check_is_mobster_mvdata = function(x)
{
  if(!inherits(x, 'mbst_data'))
  {
    stop('TODO Error message - oggetto non di classe')

  }
}
