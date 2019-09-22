#' Summary for mvMOBSTER objects.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ...
#'
#' @return Nothing
#'
#' @export
#'
#' @examples
#' TODO
summary.mbst_data = function(x, ...)
{
  print_header_data(x)

  pioTit('Mutation data')

  pioStr("\nDepth (DP)", '', suffix = '\n')

  DP(x) %>%
    group_by(sample, variable) %>%
    summarize(
      mean = round(mean(value), 2),
      median = round(median(value), 2),
      min = round(min(value), 2),
      max = round(max(value), 2)
    ) %>%
    ungroup %>%
    pioDisp

  pioStr("\nNumber of reads with variant (NV)", '', suffix = '\n')

  NV(x) %>%
    group_by(sample, variable) %>%
    summarize(
      mean = round(mean(value), 2),
      median = round(median(value), 2),
      min = round(min(value), 2),
      max = round(max(value), 2)
    ) %>%
    ungroup %>%
    pioDisp

  pioStr("\nVariant allele frequency (VAF)", '', suffix = '\n')

  VAF(x) %>%
    group_by(sample, variable) %>%
    summarize(
      mean = round(mean(value), 2),
      median = round(median(value), 2),
      min = round(min(value), 2),
      max = round(max(value), 2)
    ) %>%
    ungroup %>%
    pioDisp

  if (!is.null(x$mutations_annotations))
    pio::pioStr('Annotations', prefix = '\n', suffix = '\n',
                paste0(setdiff(colnames(x$mutations_annotations), 'id'),
                       collapse = ', '))

  # CNA data
  pioTit('Copy Number Alteration data')
  print(x$CNAqc)

  #

  print('stats about the fits...?')



}

#' Print for mvMOBSTER objects.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ...
#'
#' @return Nothing
#'
#' @export
#'
#' @examples
#' TODO
print.mbst_data = function(x, ...)
{
  print_header_data(x)

  # Print MOBSTER fits
  if (has_mobster_fits(x)) {
    lapply(names(x$fit_MOBSTER), function(w) {
      pioTit('MOBSTER fit for sample', w)
      print(x$fit_MOBSTER[[w]]$best)
    })
  }

  # Print VIBER fits
  if (has_viber_fits(x)) {
    print(x$fit_VIBER)
  }

  # Logged operations
  pio::pioTit('LOGGED OPERATIONS')

  print(x$operationLog)
}

print_header_data = function(x)
{
  check_is_mobster_mvdata(x)

  pioHdr(
    paste("mvMOBSTER dataset"),
    toPrint =
      c(
        `  Dataset` = x$description,
        `  Samples` = paste(x$samples, collapse = ', '),
        ` Purities` = paste(x$purity, collapse = ', '),
        `Mutations` = paste('N =', N(x))
      )
  )

  cat('\n')
  pioStr("MOBSTER analysis", has_mobster_fits(x), suffix = '\n')
  pioStr("  VIBER analysis", has_viber_fits(x), suffix = '\n')

}

