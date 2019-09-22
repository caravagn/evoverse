#' Filter all mutations in a VAF range
#'
#' @description Filter all mutations that, when they have
#' VAF above 0, they are not in a certain range, given by
#' \code{min.vaf}, \code{max.vaf}.
#
#' @param x A dataset of class \code{"mbst_data"}.
#' @param min.vaf Minimum VAF.
#' @param max.vaf Maximum VAF
#'
#' @return A dataset with no mutations outside the
#' required VAF range.
#' @export
#'
#' @examples
#' TODO
filter_vafrange = function(x, min.vaf, max.vaf)
{
  check_is_mobster_mvdata(x)

  pioStr("Remove mutations with VAF outside range ", min.vaf, ' - ', max.vaf)

  ids = VAF(x) %>%
    filter(value > 0 &
             (value > min.vaf | value < max.vaf)) %>%
    select(id) %>%
    distinct() %>%
    pull(id)

  all_ids = keys(x)

  pioStr("Mutations to remove ", length(setdiff(all_ids, ids)))

  x = delete_entries(x, setdiff(all_ids, ids))

  # Log update
  x = logOp(x, paste0("Entries with VAF outside range ", min.vaf, '-', max.vaf, ' removed.'))

  x
}
