#' Create a dataset.
#'
#' @description Create a mvMOBSTER dataset that can contain
#' mutation data for a one, or more, bulk samples of the same
#' tumour. The method requires mutations, copy number segments
#' and sample purity to adjust the raw VAF for copy state and
#' tumour content. The adjusted VAF of clonal mutations in
#' diploid region is expect to be around 0.5 (heterozygous
#' mutation) that corresponds to 100% of cancer cells fraction.
#' Basic filtering options are provided.
#'
#' @param mutations Dataframe of mutations. Must contain mutation
#' coordinate as chr, from and to. Plus, for every sample named
#' \code{x}, it must contain columns \code{x.DP}, \code{x.NV}
#' and \code{x.VAF} with the values of tumour depth (DP), number
#' of reads with the mutant allele (NV) and variant allele
#' frequency (VAF). Exra columns are retained.
#' @param samples Vector of samples names.
#' @param segments Copy number segments, must contain segment
#' coordinate as chr, from and to. for every sample named
#' \code{x}, it must contain columns \code{x.minor} and
#' \code{x.Major} with the minor and major number of copies
#' of the segment.
#' @param purity Vector of purities, with samples as names.
#' @param description Dataset synopsis.
#'
#' @return An object that represents a dataset from class
#' \code{mbs_data}.
#'
#' @import CNAqc
#'
#' @export
#'
#' @examples
#' data(example_mvmobster)
#' mutations = example_mvmobster$mutations
#' segments = example_mvmobster$segments
#' purity = example_mvmobster$purity
dataset = function(
  mutations,
  segments,
  samples,
  purity,
  description = "My multi-region MOBSTER dataset"
)
{
  # Structures to create
  # - mutations: query-ready tibble with VAF, DP and NV, indexable by sample, with id key
  # - mutations: location map, with id key
  # - mutations: all-remaining information(s), indexable by sample, with id key
  # - CNA: CNAqc objects

  pio::pioHdr("mvMOBSTER dataset")

  # Convert input if it is in the form of list
  if(is.list(mutations))
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
