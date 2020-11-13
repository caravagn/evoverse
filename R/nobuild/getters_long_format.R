# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plain getters (don't transform the data) in long-format
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#' Get the keys used to index the mutations in the data.
#'
#' @description Every mutation has its own key which is used to
#' access and track the mutation inside an  \code{evoverse} object.
#'
#' With this function you can extract the keys:  many functions accept
#' them as parameter to extract data.
#'
#' @param x An `evoverse` object.
#'
#' @return The keys that index stored mutations.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' my_keys = keys(example_evoverse)
#'
#' # Example usage
#' VAF(example_evoverse, ids = my_keys[1])
keys = function(x) {
  unique(x$mutations$id)
}

#' Extract VAF values as they are stored inside an \code{evoverse} dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the VAF values. These can be subset by sample and mutation
#' id; by default all entries are returned. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Default is everything
#' VAF(example_evoverse)
#'
#' # Only one sample
#' VAF(example_evoverse, samples = example_evoverse$samples[1])
VAF = function(x,
               ids = keys(x),
               samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'VAF' &
             id %in% ids &
             sample %in% samples)
}

#' Extract depth values as they are stored inside an \code{evoverse} dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the DP values. These can be subset by sample and mutation
#' id; by default all entries are returned. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Default is everything
#' DP(example_evoverse)
#'
#' # Only one sample
#' DP(example_evoverse, samples = example_evoverse$samples[1])
DP = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'DP' &
             id %in% ids &
             sample %in% samples)
}

#' Extract the number of reads with the alternative allele, as they are stored inside an \code{evoverse} dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the NV values. These can be subset by sample and mutation
#' id; by default all entries are returned. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x An `evoverse` object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries.
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Default is everything
#' NV(example_evoverse)
#'
#' # Only one sample
#' NV(example_evoverse, samples = example_evoverse$samples[1])
NV = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'NV' &
             id %in% ids &
             sample %in% samples)
}
