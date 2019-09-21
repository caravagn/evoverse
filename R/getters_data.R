#' Get the keys to access mutations.
#'
#' @description Every mutation has its own key which is used to
#' access and track the mutation. With this function you can
#' extract the keys. Many functions accept them as parameter to
#' extract data.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#'
#' @return The keys that id mutations.
#'
#' @export
#'
#' @examples
#' TODO
keys = function(x) {
  unique(x$mutations$id)
}

#' Extract VAF values from a mvMOBSTER dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the VAF values. These can be subset by sample and mutation
#' id; by default all entries are returned. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The id of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
#' TODO
VAF = function(x,
               ids = keys(x),
               samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'VAF' &
             id %in% ids &
             sample %in% samples)
}

#' Extract depth values from a MOBSTER dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the DP values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#'
#' @export
#'
#' @examples
#' TODO
DP = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'DP' &
             id %in% ids &
             sample %in% samples)
}

#' Extract the number of reads with the alternative allele from a MOBSTER dataset.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the NV values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
#' TODO
NV = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$mutations %>%
    filter(variable == 'NV' &
             id %in% ids &
             sample %in% samples)
}

#' Extract annotation values.
#'
#' @description
#'
#' Extract from the internal representation of an object all the entries
#' that refer to the annotation values. These are all values which are not
#' any of VAF, DP or NV, or any of the location coordinates for the input
#' mutations. As for the other getters, you can subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
#' TODO
Annotations = function(x,
                       ids = keys(x))
{
  x$mutations_annotations %>%
    filter(id %in% ids)
}

#' Tabular VAF values.
#'
#' @description
#'
#' Similarly to \code{VAF}, this function however returns a spread-like table
#' with one column per VAF value. As \code{VAF} the entries can be subset as
#' required, and the columns named appending a custom suffix to sample names.
#'
#' @param x A MOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.VAF} when \code{suffix = '.VAF'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' TODO
VAF_table = function(x,
                     ids = keys(x),
                     samples = x$samples,
                     suffix = '.VAF')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = VAF(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, VAF)
    output = full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}

#' Tabular NV values.
#'
#' @description
#'
#' Similarly to \code{NV}, this function however returns a spread-like table
#' with one column per VAF value. As \code{NV} the entries can be subset as
#' required, and the columns named appending a custom suffix to sample names.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.NV} when \code{suffix = '.NV'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' TODO
NV_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    suffix = '.NV')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = NV(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, NV)

    output = full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}


#' Tabular DP values.
#'
#' @description
#'
#' Similarly to \code{DP}, this function however returns a spread-like table
#' with one column per VAF value. As \code{DP} the entries can be subset as
#' required, and the columns named appending a custom suffix to sample names.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.DP} when \code{suffix = '.DP'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' TODO
DP_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    suffix = '.DP')
{
  output = tibble(`id` = ids)

  for(s in samples) {
    entry = DP(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, DP)

    output = full_join(output, entry, by = 'id')
  }

  colnames(output) = c('id', paste0(samples, suffix))

  output
}


#' Tabular VAF, DP and NV values, with annotations.
#'
#' @description
#'
#' This is just a wrapper to a combined call of function \code{VAF_table}, \code{DP_table} and
#' \code{NV_table}. The resulting outputs are bound by column; ususal subset options are available.
#' The output can be augmented with annotations from each available mutation.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param annotations False by default; if true annotations are also returned.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
#' TODO
Data_table = function(x,
                      ids = keys(x),
                      samples = x$samples,
                      annotations = FALSE)
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

#' Return the number of mutations in the dataset.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#'
#' @return The number of mutations in the dataset.
#' @export
#'
#' @examples
#' TODO
N = function(x) {
  nrow(VAF(x, samples = x$samples[1]))
}

#' Return the number of samples in the dataset.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#'
#' @return The number of samples in the dataset.
#' @export
#'
#' @examples
#' TODO
S = function(x) {
  length(x$samples)
}


