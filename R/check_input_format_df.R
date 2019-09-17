# Check a data.frame input
check_input_format_df = function(mutations,
                              segments,
                              samples,
                              purity)
{
  if(!is.data.frame(mutations)) stop("Mutations must be in a dataframe")
  if(!is.data.frame(segments)) stop("Segments must be in a dataframe")

  if(!is.vector(samples)) stop("Samples names must be in a vector")
  if(!is.vector(purity)) stop("Purity values must be in a vector")

  if(!all(samples %in% names(purity))) stop("Missing purity value for some samples")

  # Mutation columns
  muts.columns = c('chr', 'from', 'to', 'ref', 'alt')
  DP.columns = paste0(samples, ".DP")
  NV.columns = paste0(samples, ".NV")
  VAF.columns = paste0(samples, ".VAF")

  mut_all.columns = c(muts.columns, DP.columns, NV.columns, VAF.columns)

  if(any(!(mut_all.columns %in% colnames(mutations))))
    stop("Missing columns from mutation data. Required: ", paste(mut_all.columns, collapse = ', '))

  # CNA columns
  segments.columns = c('chr', 'from', 'to')
  Major.columns = paste0(samples, ".Major")
  minor.columns = paste0(samples, ".minor")

  seg_all.columns = c(segments.columns, Major.columns, minor.columns)

  if(any(!(seg_all.columns %in% colnames(segments))))
    stop("Missing columns from CNA data. Required: ", paste(seg_all.columns, collapse = ', '))

  return(
    list(
      mut_all.columns,
      seg_all.columns
    )
  )
}
