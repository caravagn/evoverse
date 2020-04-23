#' Filter mutations with VAF outside a range.
#'
#' @description This function filters mutations by VAF.
#'
#' The filter deletes all mutations that, when they have
#' VAF above 0 in a sample (and so they are "called" in that sample),
#' have a VAF value outside a desired range.
#'
#' @param x An `evoverse` object.
#' @param min.vaf Minimum VAF.
#' @param max.vaf Maximum VAF.
#'
#' @return A dataset with no mutations outside the
#' required VAF range.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' N(example_evoverse)
#'
#' N(filter_vafrange(example_evoverse, .1, .6))
filter_vaf_range = function(x, min.vaf, max.vaf)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with VAF outside range ", min.vaf, ' - ', max.vaf, suffix = '\n')

  # Group by mutation, take entries > 0 and compute min and max
  ids_to_cancel = VAF(x) %>%
    group_by(id) %>%
    filter(value > 0) %>%
    summarise(m = min(value), M = max(value)) %>%
    filter(m < min.vaf | M > max.vaf) %>%
    pull(id) %>%
    unique

  pioStr("Mutations to remove ", length(ids_to_cancel), suffix = '\n')

  x = delete_entries(x, ids_to_cancel)

  # Log update
  x = logOp(x, paste0("Entries with VAF outside range ", min.vaf, '-', max.vaf, ' removed.'))

  x
}
