




# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Private auxiliary functions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)

  obj
}



# Conditional mutate
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}



