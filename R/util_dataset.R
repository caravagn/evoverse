# Log an operation
logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)

  obj
}

# Check if the object has MOBSTER fits
has_mobster_fits = function(x)
{
  !all(is.null(x$fit_MOBSTER))
}

# Check if the object has VIBER fits
has_viber_fits = function(x)
{
  !all(is.null(x$fit_VIBER))
}

# Check if the object is of class evoverse, raise error if not
check_is_mobster_mvdata = function(x)
{
  if(!inherits(x, 'mbst_data'))
  {
    stop(
      'The input object is not a evoverse dataset, aborting.\nHead\n',
      head(str(x))
      )
  }
}
