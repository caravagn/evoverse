#' Data analysis template
#'
#' @description
#'
#' A three-steps higher-order pipeline; loads input data, computes a result and
#' postprocess it. It implements an cache system to avoid computing when results
#' are available (optional).
#'
#' In tidy language the pipeline is roughly this, plus the caching system
#'
#' fun_data_load %>%
#'    fun_compute %>%
#'    fun_postprocess
#'
#' @param x A sample character id.
#' @param fun_data_load A function applied to sample to obtain an input object;
#' for instance loads a CSV from a folder indexed by the sample id.
#' @param fun_compute A function that implements the actual computation on
#' to the input object; for instance it clusters the input. All parameters must
#' be made explicit inside this function. .
#' @param fun_postprocess A function to compute the final result from the
#' analysis output, and its input.
#' @param cache By default \code{NULL}, which turns off the cache. If this is a
#' filename, and the file exists, no computation is skip.
#'
#' @return A list with the sample and results from each each step of the
#' pipeline
#' @export
#'
#' @examples
evo_pipe = function(x,
                    fun_data_load,
                    fun_compute,
                    fun_postprocess,
                    cache = NULL)
{
  stopifnot(is.character(x))
  stopifnot(is.function(fun_data_load))
  stopifnot(is.function(fun_compute))
  stopifnot(is.function(fun_postprocess))

  # Stopwatch
  TIME_all = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  difftime_ = function(TIME)
  {
    TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
                    TIME, units = "mins")
    return(format(TIME, digits = 0))
  }

  # Startup
  cli::cli_rule("{.field {x}}",
                right = format(Sys.time(), "%a %b %d %X %Y"))

  m = function(x)
    crayon::bgGreen(crayon::black(x))
  o = function(x)
    crayon::bgMagenta(crayon::white(x))

  # 1 - Input

  # Cache output check
  if (!is.null(cache) && file.exists(cache))
  {
    cli::cli_alert_warning("{o(' Cache output found ')} will not compute the analysis. ({.field {x}})")
    stop("Cached!")
  }

  # Load task input
  TIME_input = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  cli::cli_alert_info("{m(' Loading input ')} {.field {x}}")

  input = fun_data_load(x)

  cli::cli_alert_success(
    "{m(' Loading input ')} completed for {.field {x}} in {.value {difftime_(TIME_input)}}"
  )

  # 2 - Computation

  TIME_compute = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  cli::cli_alert_info("{m(' Computing analysis ')} {.field {x}}")

  analysis = fun_compute(input)

  cli::cli_alert_success(
    "{m(' Computing analysis ')} completed for {.field {x}} in {.value {difftime_(TIME_compute)}}"
  )

  # 3 - Post-process results

  TIME_postprocess = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  cli::cli_alert_info("{m(' Processing results ')} {.field {x}}")

  postprocess = fun_postprocess(input, analysis)

  cli::cli_alert_success(
    "{m(' Processing results ')} completed for {.field {x}} in {.value {difftime_(TIME_postprocess)}}"
  )

  # Bye
  cli::cli_rule(
    "{crayon::bold('eoverse')} completed {.field {x}} in {.value {difftime_(TIME_all)}}",
    right = format(Sys.time(), "%a %b %d %X %Y")
  )

  return(list(
    x = x,
    input = input,
    analysis = analysis,
    postprocess = postprocess
  ))
}
