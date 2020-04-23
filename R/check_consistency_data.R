check_consistency_data = function(x)
{
  list_mutations = keys(x)
  dt = Data_table(x)

  # No columns with numeric infos is NA
  cannot_be_NA = dt %>%
    select(ends_with('VAF'), ends_with('DP'), ends_with('NV')) %>%
    unlist %>%
    is.na %>%
    any

  if(cannot_be_NA)
    message("There are NA values into some entries for VAF, DP and NV -- there should not be such values.")

  # No mutation DP is 0
  cannot_be_0_DP = dt %>%
    select(ends_with('DP')) %>%
    unlist
  cannot_be_0_DP = any(cannot_be_0_DP == 0)

  if(cannot_be_0_DP)
    message("There are 0 values into some entries for DP -- there should not be such values (coverage should be reported even if the mutation is not called in a sample).")

  return(
    any(cannot_be_NA, cannot_be_0_DP)
  )
}


