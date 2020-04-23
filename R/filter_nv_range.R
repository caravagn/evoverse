#' Filter mutations outside a range of NV values.
#'
#'  @description This function filters mutations by NV (number of variants).
#'
#' The filter deletes all mutations that, at least in one sample,
#' a number of variants NV outside a desired range.
#'
#' @param x An `evoverse` object.
#' @param min.NV Minimum NV
#' @param max.NV Maximum NV
#'
#' @return A dataset with no mutations outside the
#' required NV range.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' N(example_evoverse)
#'
#' N(filter_nv_range(example_evoverse, 50, 100))
filter_nv_range = function(x, min.NV, max.NV)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with NV outside range ", min.NV, ' - ', max.NV, suffix = '\n')

  # Group by mutation, take entries > 0 and compute min and max
  ids_to_cancel = NV(x) %>%
    group_by(id) %>%
    filter(value > 0) %>%
    summarise(m = min(value), M = max(value)) %>%
    filter(m < min.NV | M > max.NV) %>%
    pull(id) %>%
    unique

  pioStr("Mutations to remove ", length(ids_to_cancel), suffix = '\n')

  x = delete_entries(x, ids_to_cancel)

  # Log update
  x = logOp(x, paste0("Entries with NV outside range ", min.NV, '-', max.NV, ' removed.'))

  x
}
