... = function(x)
{



  # Take segments and count how many mutations they have mapped
  pass_segments = lapply(x$CNAqc,
                         function(x)
                           x$cna %>%
                           mutate(n = ifelse(is.na(n), FALSE, n > N.min)) %>%
                           filter(n) %>%
                           pull(segment_id)
  )
  pass_segments = Reduce(intersect, pass_segments)

  # A summary table with the number of mutations per segment
  size_table = x$map_mut_seg %>%
    group_by(seg_id) %>%
    summarise(N = length(unique(id)))

  # we reject those with N < N.min
  rejected = size_table %>%
    filter(N < !!N.min) %>%
    arrange(desc(N)) %>%
    inner_join(y = x$segments, by=c("seg_id" = "id")) %>%
    select(seg_id, N, chr, from, to) %>%
    distinct

  accepted = size_table %>%
    filter(N >= !!N.min) %>%
    arrange(desc(N)) %>%
    inner_join(y = x$segments, by=c("seg_id" = "id")) %>%
    select(seg_id, N, chr, from, to) %>%
    distinct

  pioTit(paste0("Segments report (cutoff " , N.min, ' muts/seg)'))

  pioStr(paste0("\nN < ", N.min), nrow(rejected), suffix = '(rejected)\n')
  print(rejected)

  pioStr(paste0("\nN >= ", N.min), nrow(accepted), suffix = '(accepted)')
  print(accepted)

  if(nrow(accepted) == 0)
    stop("No segments can be used, aborting.")

  # Subset to match accepted ... Non-accepted mutations can be immediately removed
  rejected = x$map_mut_seg %>%
    filter(seg_id %in% rejected$seg_id) %>%
    pull(id)

  # N(x)

  if(length(rejected) > 0)
    x = delete_entries(x, unique(rejected))

  pio::pioTit(paste0("Adjusting the VAF for ",
                     nrow(x$map_mut_seg)," entries across ", length(x$samples), " samples"))

  new.VAF = NULL

  # for every segment to map
  seg_to_match = unique(x$map_mut_seg$seg_id)

  if(nrow(segments) > 1) pb = txtProgressBar(min = 0, max = length(seg_to_match), style = 3)

  for(seg in seq(seg_to_match))
  {
    if(length(seg_to_match) > 1) setTxtProgressBar(pb, seg)

    # for every sample
    for(s in x$samples)
    {
      # get VAF values for these entries
      values = VAF(x,
                   ids = x$map_mut_seg %>% filter(seg_id == seg_to_match[seg]) %>% pull(id),
                   samples = s
      ) %>% # Adjust VAF according to minor/ Major (CN=m+M)
        mutate(
          Major = Major(x, seg_to_match[seg], s),
          minor = minor(x, seg_to_match[seg], s),
          CN = Major + minor,
          purity = purity[s],
          adj_VAF = vaf_adjustCN(value, minor, Major, purity)
        )

      new.VAF = bind_rows(new.VAF, values)
    }
  }

  # Save a copy of the mapping
  x$VAF_cn_adjustment = new.VAF

  new.VAF$value = new.VAF$adj_VAF
  new.VAF = new.VAF %>% select(id, variable, value, sample)

  x$data = bind_rows(
    new.VAF,
    DP(x),
    NV(x)
  )

  # As CN, we can keep only those accepted
  x$segments = x$segments %>%
    filter(id %in% accepted$seg_id)

  # Log update
  x = logOp(x, "Added CN and adjusted VAF")

  all_zeroes(x)
}
