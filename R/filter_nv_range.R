#' Filter all mutations in a NV range
#'
#' @description Filter all mutations that, when they have
#' DP above 0, they are not in a certain range, given by
#' \code{min.NV}, \code{max.NV}.
#
#' @param x A dataset of class \code{"mbst_data"}.
#' @param min.NV Minimum NV
#' @param max.NV Maximum NV
#'
#' @return A dataset with no mutations outside the
#' required NV range.
#' @export
#'
#' @examples
#' TODO
filter_dp_range = function(x, min.NV, max.NV)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with NV outside range ", min.NV, ' - ', max.NV)

  ids = NV(x) %>%
    filter(value > 0 &
             (value > min.NV | value < max.NV)) %>%
    select(id) %>%
    distinct() %>%
    pull(id)

  all_ids = keys(x)

  pioStr("Mutations to remove ", length(setdiff(all_ids, ids)))

  x = delete_entries(x, setdiff(all_ids, ids))

  # Log update
  x = logOp(x, paste0("Entries with NV outside range ", min.NV, ' - ', max.NV, ' removed.'))

  x
}
