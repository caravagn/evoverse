#' Filter mutations outside a range of coverage values.
#'
#'  @description This function filters mutations by DP (depth).
#'
#' The filter deletes all mutations that, at least in one sample,
#' have read depth DP outside a desired range.
#'
#' @param x An `evoverse` object.
#' @param min.DP Minimum DP
#' @param max.DP Maximum DP
#'
#' @return A dataset with no mutations outside the
#' required DP range.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' N(example_evoverse)
#'
#' N(filter_dp_range(example_evoverse, 50, 200))
filter_dp_range = function(x, min.DP, max.DP)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with DP outside range ", min.DP, ' - ', max.DP, suffix = '\n')

  # Group by mutation, take entries > 0 and compute min and max
  ids_to_cancel = DP(x) %>%
    group_by(id) %>%
    filter(value > 0) %>%
    summarise(m = min(value), M = max(value)) %>%
    filter(m < min.DP | M > max.DP) %>%
    pull(id) %>%
    unique

  pioStr("Mutations to remove ", length(ids_to_cancel), suffix = '\n')

  x = delete_entries(x, ids_to_cancel)

  # Log update
  x = logOp(x, paste0("Entries with DP outside range ", min.DP, '-', max.DP, ' removed.'))

  x
}
