# Special getter

#' Extract all the annotations available for the mutations stored inside an \code{evoverse} dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that do not refer to any of VAF, DP and NV values, chromsome locations.
#'
#' This getter differs from the one for DP, NV and VAF values. In particular
#' since annotations are the same across all samples, this function does not
#' have samples as parameters. Moreover, in this case the data is stored
#' in the wide format inside \code{evoverse}, and this functions returns
#' that type of representation.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#'
#' @return A tibble of the required entries.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Default is everything
#' Annotations(example_evoverse)
#'
#' # Only one sample
#' Annotations(example_evoverse)
Annotations = function(x,
                       ids = keys(x))
{
  x$mutations_annotations %>%
    filter(id %in% ids)
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Wide getters (transform the data) for wide-format
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#' Extract VAF values in wide format (multi-column).
#'
#' @description
#'
#' This function spreads the result of the call to the \link{\code{VAF}}
#' getter following the principles of wide data representation.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.VAF} when \code{suffix = '.VAF'}.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' VAF_table(example_evoverse)
VAF_table = function(x,
                     ids = keys(x),
                     samples = x$samples,
                     suffix = '.VAF')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = VAF(x, ids = ids, samples = s) %>% spread(variable, value) %>% dplyr::select(id, VAF)
    output = dplyr::full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}

#' Extract NV values in wide format (multi-column).
#'
#' @description
#'
#' This function spreads the result of the call to the \link{\code{NV}}
#' getter following the principles of wide data representation.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.VAF} when \code{suffix = '.NV'}.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' NV_table(example_evoverse)
NV_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    suffix = '.NV')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = NV(x, ids = ids, samples = s) %>% spread(variable, value) %>% dplyr::select(id, NV)

    output = dplyr::full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}


#' Extract DP values in wide format (multi-column).
#'
#' @description
#'
#' This function spreads the result of the call to the \link{\code{DP}}
#' getter following the principles of wide data representation.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.VAF} when \code{suffix = '.DP'}.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' DP_table(example_evoverse)
DP_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    suffix = '.DP')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = DP(x, ids = ids, samples = s) %>% spread(variable, value) %>% dplyr::select(id, DP)

    output = full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}


#' Extract all the data stored in an \code{evoverse} dataset, in the wide format.
#
#' @description
#'
#' This is just a wrapper to a combined call of function \link{\code{VAF_table}}, \link{\code{DP_table}} and
#' \link{\code{NV_table}}. The resulting outputs are joind by id; usual subsetting options are available.
#' The output contains also the annotations from each available mutation, obtained with \link{\code{Annotations}}.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' Data_table(example_evoverse)
Data_table = function(x,
                      ids = keys(x),
                      samples = x$samples)
{
  # Joined data tables
  dt = full_join(
    full_join(
      VAF_table(x, ids, samples),
      DP_table(x, ids, samples),
      by = 'id'),
    NV_table(x, ids, samples),
    by = 'id'
  )

  # All variables
  an = Annotations(x, ids)

  full_join(dt, an, by = 'id')
}
