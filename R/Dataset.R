#' Create an \code{evoverse} dataset.
#'
#' @description An \code{evoverse} dataset contains, for one or more bulk samples
#' of the same tumour, the following sequencing information:
#' \enumerate{
#' \item mutation data, in the usual format of substitutions involving
#' one or more nucleotides;
#' \item absolute copy number data, reported as segments profiles;
#' \item sample purity (proportion of tumour cells in each one of the samples).
#' }
#' This method uses all these information to create an object of
#' class \code{mbs_data}, which provides different S3 methods to
#' print, visualize and query the loaded data.
#'
#' Two different types of data formats for mutation and copy number data
#' are supported by this function, which will try to determine
#' automatically the type of input. All other parameters have a trivial
#' data type, independently of mutation and copy number segments data.
#'
#' Format one for mutations and copy number segments (minimum information
#' required):
#'
#' \itemize{
#' \item A `dataframe` that must contain mutation coordinates as `chr`,
#' `from`, `to` in `hg19` chromosome coordinates, plus `ref` and `alt`
#' for the reference and alternative alleles. For every sample `X`,
#' this dataframe must contain columns `X.DP`, `X.NV` and `X.VAF` for the depth
#' (total number of reads with both alleles), the number of reads
#' with the alternative allele and the allele frequency (\code{VAF = NV/DP}).
#' Chromosome names must be in the format `chr1`, `chr2`, etc.; alleles should
#' be characters and all other fields numeric.
#' \item A `dataframe` that must contain copy number segments as `chr`,
#' `from`, `to` in `hg19` chromosome coordinates, plus `X.Major` and `X.minor`
#' for the number of copies of the major allele, and the minor (B-allele)
#' allele in sample `X`'s segment. Chromosome names must be in the same
#' format used by mutations.
#' }
#'
#' Format two for mutations and copy number segments (minimum information
#' required):
#'
#' \itemize{
#' \item  A list of `dataframes` that contains mutations for a unique sample
#' (instead than for all of the samples). Mutation genomic coordinates
#' must follow the same guidelines and format discussed above, but this time
#' the list entry for patient `X` can have columns `DP`, `NV` and `VAF`
#' with the usual meaning. The list should be named, therefore the dataframe
#' for sample `X` should be named `X`.
#' \item  A list of `dataframes` that contains segments. As for mutations, in the list
#' named for patient `X` the major an minor alleles can just be named `Major` and
#' `minor` (i.e., without the sample appearing the column name).
#' }
#'
#' Regardless the format, mutations are identified using as key the chromosome
#' name, the genome location, the reference and alternative alleles. For instance,
#' a SNV with a C>T substitution at chromosome 6, and nucleotide 10048079 will have
#' id \code{"chr6:10048079:10048079:C:T"}. The ids are used to create an internal
#' representation of the data that allows easy querys and modifications; and are
#' accessibly with function \link{\code{keys}}. Extra columns
#' available in the mutation and copy number data will be retained in internal
#' representation of the data, and won't be manipulated.
#'
#' \strong{Important consistency requirements:} this function makes basic checks of
#' input consistency, flagging dubious inputs and in some cases refusing to proceed.
#' It is important that the user checks that:
#' \itemize{
#' \item All mutations are annotated in every sample, even if absent. So even if a mutation
#' `x``has VAF 0 in a sample, its actual depth, NV, and VAF should be annotated.
#' \item Copy number segments should be spanning the same intervals across all the input
#' samples. So this means that segments calls are joinly called even if the samples are
#' analysed independently.
#' }
#'
#'\strong{Important note:} this package uses the \code{CNAqc} package to map mutations
#' to copy number segments, and processing/ visualising copy number calls. Please see the
#' webpage \url{https://github.com/caravagn/CNAqc} for how to use this
#' package to asses the quality of your copy number calls.
#'
#' @param mutations Dataframe of mutations in all samples, or list of dataframes of
#' mutations, one entry per sample. The formats and requirements are describe above.
#' @param segments Dataframe of copy number segments in all samples, or list of dataframes of
#' copy number segments, one entry per sample. The formats and requirements are describe above.
#' @param samples Vector of samples names.
#' @param purity Vector of purities, named consistently to \code{samples}.
#' @param description Dataset synopsis.
#'
#' @return An object that represents a dataset from class \code{mbs_data}.
#'
#' @seealso	Two related packages developed by Giulio Caravagna
#' \itemize{
#' \item MOBSTER, the model-based approach to cluster tumour sequencing
#' data, available at \url{https://github.com/caravagn/MOBSTER}.
#' \item CNAqc, the Copy Number Analysis quality check package,
#'  available at \url{https://github.com/caravagn/CNAqc}.
#' }
#'
#' @import CNAqc
#' @import mobster
#' @import pio
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#' mutations = example_mvmobster$mutations
#' segments = example_mvmobster$segments
#' purity = example_mvmobster$purity
dataset = function(
  mutations,
  segments,
  samples,
  purity,
  description = "My evoverse dataset"
)
{
  # Structures to create
  # - mutations: query-ready tibble with VAF, DP and NV, indexable by sample, with id key
  # - mutations: location map, with id key
  # - mutations: all-remaining information(s), indexable by sample, with id key
  # - CNA: CNAqc objects

  pio::pioHdr("evoverse dataset")

  # Convert input if it is in the form of list (not tibble)
  if(is.list(mutations) & !is_tibble(mutations))
  {
    conversion =  converted_dataset(
      mutations,
      segments,
      purity,
      samples
    )

    mutations = conversion[[1]]
    segments = conversion[[2]]
    purity = conversion[[3]]
    samples = conversion[[4]]
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check input format for data
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  columns_required =
    check_input_format_df(mutations,
                          segments,
                          samples,
                          purity)

  mutation_columns_required = columns_required[[1]]
  segments_columns_required = columns_required[[2]]

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Create IDs for database, and split information to retain in different formats
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  mutations = mutations %>%
    mutate(id = paste(chr, from, to, ref, alt, sep = ':'))

  # retain only the information we actually need
  mutations_retained = mutations %>% select(-!!mutation_columns_required, id)
  mutations = mutations %>% select(!!mutation_columns_required, id)
  mutations_locations = mutations %>% select(chr, from, to, ref, alt, id)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # NA values in the data are notified
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # NAs in any of the columns for values
  if(any(is.na(mutations)))
  {
    warning("NA values in the DP, NV and VAF fields. Some mutations will be removed, you might want to check your data....")

    mutations = mutations[
      complete.cases(mutations), , drop = FALSE
      ]
  }

  pioStr('Mutations', paste0('N = ', nrow(mutations)), suffix = '\n')

  require(CNAqc)

  # Sample-specific CNA maps data
  CNAqc_mappings =
    lapply(
      samples,
      function(s)
      {
        cat('\n')
        pioTit(
          "Creating CNAqc object for", s
        )

        data_s = mutations %>% select(chr, from, to, ref, alt, id, starts_with(paste0(s, '.')))
        colnames(data_s) = gsub(paste0(s, '\\.'), '', colnames(data_s))
        colnames(data_s) = gsub('\\.', '', colnames(data_s))

        segments_s = segments %>% select(chr, from, to, starts_with(paste0(s, '.')))
        colnames(segments_s) = gsub(paste0(s, '\\.'), '', colnames(segments_s))
        colnames(segments_s) = gsub('\\.', '', colnames(segments_s))

        CNAqc::init(data_s, segments_s, purity[s])
  })
  names(CNAqc_mappings) = samples

  # Mapped mutations from CNAqc
  mapped_mutations = lapply(
    names(CNAqc_mappings),
    function(x)
      CNAqc_mappings[[x]]$snvs %>%
      select(VAF, DP, NV, karyotype, id) %>%
      mutate(sample = x) %>%
      reshape2::melt(id = c('karyotype', 'id', 'sample')) %>%
      mutate(variable = paste(variable)) %>%
      as_tibble
    )

  mapped_mutations = Reduce(bind_rows, mapped_mutations) %>%
    select(id, sample, variable, value, karyotype)

  # Non-mappable mutation are cleaned up (removed across all samples)
  non_mappable = mapped_mutations %>%
    group_by(id, sample) %>%
    filter(is.na(karyotype)) %>%
    pull(id) %>%
    unique

  # notify which wkll be removed
  cat('\n')
  pioTit('Non-mappable mutations')
  mapped_mutations %>%
    filter(id %in% non_mappable, variable == 'VAF') %>%
    spread(variable, value) %>%
    spread(sample, id) %>%
    pioDisp

  # Actual clean up of the data
  mutations = mapped_mutations %>%
    filter(!(id %in% non_mappable))

  mutations_locations = mutations_locations %>%
    filter(!(id %in% non_mappable))

  mutations_retained = mutations_retained %>%
    filter(!(id %in% non_mappable))

  stopifnot(nrow(mutations_retained) == nrow(mutations_locations))

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Create a mbst_data object
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  x = list()
  class(x) <- "mbst_data"

  x$mutations = mutations
  x$mutations_locations = mutations_locations
  x$mutations_annotations = mutations_retained

  x$samples = samples
  x$purity = purity

  x$CNAqc = CNAqc_mappings

  x$description = description

  # Log creation
  x = logOp(x, "Initialization")

  return(x)
}
