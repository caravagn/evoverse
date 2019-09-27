#' Compute projected mutations from MOBSTER fits.
#'
#' @description After MOBSTER has been run on each sample stored in
#' \code{x}, the fits can be used to projected mutations that are
#' found to be part of a tail at least onece. For these mutations
#' their VAF and NV levels are forced to 0 in all samples, and
#' the mutations are then removed.
#'
#' @param x A dataset of class \code{"mbst_data"}.
#'
#' @return A dataset with projected mutations from MOBSTER fits.
#' @export
#'
#' @examples
#' TODO
filter_tails = function(x)
{
  check_is_mobster_mvdata(x)

  if (!has_mobster_fits(x))
  {
    warning("There are no MOBSTER clusters coputed for this dataset, will return just the object.")
    return(x)
  }

  # get clusters
  clusters = Clusters(x) %>%
    select(ends_with('.MOBSTER_cluster'), id)

  colnames(clusters) = gsub('.MOBSTER_cluster', '', colnames(clusters))

  clusters = clusters %>%
    mutate(
      is_tail = any(!!x$samples == 'Tail')
    )

  pioStr("MOBSTER clusters counts", suffix = '\n')
  print(sapply(clusters[1:(ncol(clusters) - 2)], table))

  pioStr("Tail mutations", sum(clusters$is_tail), suffix = '\n')

  to_cancel = clusters %>% filter(is_tail) %>% pull(id)
  if(length(to_cancel) > 0)
  {
    x = delete_entries(x, to_cancel)
    x = logOp(x, paste0("Removed ", sum(clusters$is_tail), " tail mutations."))
  }
  else
    message("No tail mutations to delete from this dataset.")

  return(x)
}
