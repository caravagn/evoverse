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

  colnames_ids = lapply(mutations, colnames)
  colnames_ids = Reduce(intersect, colnames_ids)
  stopifnot(all(colnames_ids %in% colnames(mutations[[1]])))

  # shared_loc_columns = c('chr', 'from', 'to', 'ref', 'alt')
  # specific_columns = c('VAF', "DP", "NV")
  # extra_columns = setdiff(colnames_ids, c(shared_loc_columns, specific_columns))

  # Add sample names to samples columns
  for(s in samples)
  {
    # Mutations
    new_vaf_name = paste0(s, '.VAF')
    new_dp_name = paste0(s, '.DP')
    new_nv_name = paste0(s, '.NV')

    mutations[[s]] = mutations[[s]] %>%
      rename(
        !!new_vaf_name := VAF,
        !!new_dp_name := DP,
        !!new_nv_name := NV
      )

    # CNA
    new_minor_name = paste0(s, '.minor')
    new_major_name = paste0(s, '.Major')

    segments[[s]] = segments[[s]] %>%
      rename(
        !!new_minor_name := minor,
        !!new_major_name := Major
      )
  }

  mutations_data = suppressMessages(
    Reduce(
      function(x, y) full_join(x, y),
      mutations)
  )

  segments_data = suppressMessages(
    Reduce(
      function(x, y) full_join(x, y),
      segments)
  )

  return(list(mutations_data, segments_data, purity, samples))
}
