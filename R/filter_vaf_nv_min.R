#' Impose minimum VAF and NV (read counts) value for mutations.
#'
#' @description This function sets the VAF and NV of mutations below a minimum cutoff to 0.
#' The mutation is not removed from the dataset.
#'
#' @param x An `evoverse` object.
#' @param min.vaf Minimum VAF.
#' @param min.NV Minimum NV
#'
#' @return A dataset where any VAF entry below \code{min.vaf}, and any NV entry below
#' \code{min.nv} are changed to 0.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' N(example_evoverse)
#'
#' N(filter_vaf_nv_min(example_evoverse, .05, 5))
#'
#' Data_table(example_evoverse)
#'
#' Data_table(filter_vaf_nv_min(example_evoverse, .05, 5))
filter_vaf_nv_min = function(x, min.vaf, min.nv)
{
  check_is_mobster_mvdata(x)

  set_min_vaf = 0
  pioStr("Changinng mutation entries to have minimum VAF/NV ", paste0(min.vaf, ' / ', min.nv, ' will become ', set_min_vaf), suffix = '\n')

  # Ids to remove
  ids_VAF = VAF(x) %>%
    filter(value > 0, value <= min.vaf) %>%
    pull(id) %>%
    unique

  ids_NV = NV(x) %>%
    filter(value > 0, value <= min.nv) %>%
    pull(id) %>%
    unique

  all_ids = c(ids_VAF, ids_NV) %>% unique
  L = length(all_ids)

  # Notification
  pioStr("Mutations to edit",  'n =', L, paste0('(', (L/N(x)) %>% round(digits = 2) , '% of ', N(x), ')'), suffix = '\n')

  # Manual modification of the entries
  x$mutations = x$mutations %>%
    mutate(
      value = ifelse((variable == 'VAF') & (value <= min.vaf) & (id %in% all_ids), set_min_vaf, value),
      value = ifelse((variable == 'NV') & (value <= min.nv) & (id %in% all_ids), set_min_vaf, value)
    )

  # Log update
  x = logOp(x, paste0("Entries with VAF/NV below ", min.vaf, ' / ', min.nv, 'set to ', set_min_vaf))

  x
}
