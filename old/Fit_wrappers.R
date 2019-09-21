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
mobster_fit_multivariate = function(x, samples = x$samples, ...)
{
  require(mobster)

  x$fit.MOBSTER = lapply(samples,
                         function(w)
                         {
                           pio::pioTit(paste("Fitting with MOBSTER sample", w))

                           data = VAF(x, samples = w) %>%
                             filter(value > 0) %>%
                             spread(variable, value)

                           print(data)

                           mf = mobster::mobster_fit(
                             x = data,
                             ...)

                           mf
                         })

  names(x$fit.MOBSTER) = samples

  cat("\n")
  pio::pioTit("MOBSTER fits")
  for(s in samples) print(x$fit.MOBSTER[[s]]$best)

  # Log update
  x = logOp(x, paste0("Fit MOBSTER to ", paste0(samples, collapse = ', ')))

  x
}


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
mobster_fit_VIBER = function(x, samples = x$samples, ...)
{
  require(VIBER)

  # Prepare inputs
  nv = NV_table(x, samples = samples)
  dp = DP_table(x, samples = samples)

  colnames(nv) = gsub(pattern = '.NV', replacement = '', colnames(nv))
  colnames(dp) = gsub(pattern = '.DP', replacement = '', colnames(dp))

  nv = nv[, colnames(dp)]

  pio::pioTit("Input for variational Binomial clustering with VIBER")

  pio::pioStr("NV values ", '')
  print(nv)
  pio::pioStr("DP values", '')
  print(dp)

  # Wrap fit call
  x$fit.Binomial = VIBER::variational_fit(
    x = nv %>% select(-id),
    y = dp %>% select(-id),
    ...
  )

  # Log update
  x = logOp(x, paste0("Fit VIBER to ", paste0(samples, collapse = ', ')))

  x
}






