#' Filter all mutations in a DP range
#'
#' @description Filter all mutations that, when they have
#' DP above 0, they are not in a certain range, given by
#' \code{min.DP}, \code{max.DP}.
#
#' @param x A dataset of class \code{"mbst_data"}.
#' @param min.DP Minimum DP
#' @param max.DP Maximum DP
#'
#' @return A dataset with no mutations outside the
#' required DP range.
#' @export
#'
#' @examples
#' TODO
filter_dp_range = function(x, min.DP, max.DP)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with DP outside range ", min.DP, ' - ', max.DP)

  ids = DP(x) %>%
    filter(value > 0 &
             (value > min.DP | value < max.DP)) %>%
    select(id) %>%
    distinct() %>%
    pull(id)

  all_ids = keys(x)

  pioStr("Mutations to remove ", length(setdiff(all_ids, ids)))

  x = delete_entries(x, setdiff(all_ids, ids))

  # Log update
  x = logOp(x, paste0("Entries with VAF outside range ", min.DP, '-', max.DP, ' removed.'))

  x
}
