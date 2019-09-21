#' Compute VIBER Binomial fits from counts data.
#'
#' @description
#'
#' This function wraps a call to the \code{variational_fit} function
#' from package \code{VIBER}, which implements a variational Binmial
#' mixture model to cluster read counts. See the package manual for
#' \code{VIBER} to see the parameters that can be forwarded, and the
#' plots that can be generated afterwards.
#'
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param samples The samples on which VIBER will be run.
#' @param ... Extra parameters that will be forwarded to
#' \code{variational_fit} function from package \code{VIBER}.
#'
#' @return A new object with a VIBER fit for the required
#' samples in \code{x$fit.Binomial}.
#'
#' @import VIBER
#'
#' @export
#'
#' @examples
#' TODO
analyze_VIBER = function(x, samples = x$samples, ...)
{
  require(VIBER)

  # Prepare inputs
  nv = NV_table(x, samples = samples)
  dp = DP_table(x, samples = samples)

  colnames(nv) = gsub(pattern = '.NV', replacement = '', colnames(nv))
  colnames(dp) = gsub(pattern = '.DP', replacement = '', colnames(dp))

  nv = nv[, colnames(dp)]

  pioTit("Variational Binomial clustering with VIBER")

  pioStr("NV values ", '', suffix = '\n')
  print(nv)
  pioStr("DP values", '', suffix = '\n')
  print(dp)

  # Wrap fit call
  x$fit_VIBER = VIBER::variational_fit(
    x = nv %>% select(-id),
    y = dp %>% select(-id),
    ...
  )

  # Store data
  x$fit_VIBER$fit_data = nv %>% select(id)

  # Log update
  x = logOp(x, paste0("Fit VIBER to ", paste0(samples, collapse = ', ')))

  x
}






