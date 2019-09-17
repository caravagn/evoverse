converted_dataset = function(
  mutations,
  segments,
  purity,
  samples
)
{
  stopifnot(all(names(mutations) == names(segments)))
  stopifnot(all(samples %in% names(mutations)))
  stopifnot(all(samples %in% names(segments)))

  # mutation_ids = lapply(
  #   mutations,
  #   function(x) {
  #     x %>%
  #       mutate(id = paste(chr, from, to, ref, alt, sep = ':'))
  # })
  # mutation_ids = Reduce(intersect, lapply(mutation_ids, function(x) x$id))

  # pio::pioStr("Indexable mutations:", length(mutation_ids))

  colnames_ids = lapply(mutations, colnames)
  colnames_ids = Reduce(intersect, colnames_ids)
  stopifnot(all(colnames_ids %in% colnames(mutations[[1]])))

  # shared_loc_columns = c('chr', 'from', 'to', 'ref', 'alt')
  # specific_columns = c('VAF', "DP", "NV")
  # extra_columns = setdiff(colnames_ids, c(shared_loc_columns, specific_columns))

  mutations_data = mutations
  segments_data = segments

  # Add sample names to samples columns
  for(s in samples)
  {
    # Mutations
    y = mutations[[s]] %>% select(chr, from, to, ref, alt, VAF, DP, NV)
    colnames(y)[6:8] = paste0(s, '.', c('VAF', 'DP', 'NV'))

    mutations_data[[s]] = y

    # CNA
    w = segments[[s]] %>% select(chr, from, to, minor, Major)
    colnames(w)[4:5] = paste0(s, '.', c('minor', 'Major'))

    segments_data[[s]] = w
  }

  mutations_data = Reduce(
    function(x, y) full_join(x, y, by = c("chr", "from", "to", "ref", "alt")),
    mutations_data)

  segments_data = Reduce(
    function(x, y) full_join(x, y, by = c("chr", "from", "to")),
    segments_data)


  return(list(mutations_data, segments_data, purity, samples))
}
