#' Compute MOBSTER fits for the available samples.
#'
#' @description
#'
#' This function wraps a call to the \code{mobster_fit} function
#' from package \code{mobster}, which implements the MOBSTER
#' statistical model. This is a combination of one Pareto power
#' law for a tumour tail, plus a mixture of Beta distributions.
#' The function applies the fitting function to the samples stored in the
#' dataset \code{x}, one by one, and the returned object contains the
#' MOBSTER fits which can be inspected with the functions available
#' in the \code{mobster} package.
#'
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param samples The samples on which MOBSTER will be run.
#' @param ... Extra parameters that will be forwarded to
#' \code{mobster_fit} function from package \code{mobster}.
#'
#' @return A new object with MOBSTER fits for the required
#' samples in \code{x$fit.MOBSTER}.
#'
#' @import mobster
#'
#' @export
#'
#' @examples
#' TODO
analyze_mobster = function(x, samples = x$samples, ...)
{
  require(mobster)

  x$fit_MOBSTER = lapply(samples,
                         function(w)
                         {
                           pio::pioTit("Fitting sample", w, "with MOBSTER")

                           V = VAF(x, samples = w)

                           if(sum(V$value == 0))
                            message(sum(V$value == 0), " mutations with 0 VAF are not fit with this sample.")

                           data = V %>%
                             filter(value > 0) %>%
                             spread(variable, value)

                           mf = mobster::mobster_fit(
                             x = data,
                             ...)

                           mf
                         })

  names(x$fit_MOBSTER) = samples

  # Log update
  x = logOp(x, paste0("Fit MOBSTER to ", paste0(samples, collapse = ', ')))

  x
}
