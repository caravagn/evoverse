byLoc = function(x, loc.id, offset_around_centromers) {
  muts_CN = x$locations %>%
    spread(variable, value)

  thisSeg = x$segments %>%
    filter(id == loc.id)

  # Enforce the correct types!
  muts_CN$chr = as.character(muts_CN$chr)
  muts_CN$from = as.numeric(muts_CN$from)
  muts_CN$to = as.numeric(muts_CN$to)

  thisSeg$chr = as.character(thisSeg$chr)
  thisSeg$from = as.numeric(thisSeg$from)
  thisSeg$to = as.numeric(thisSeg$to)

  # Within segment
  which_muts = muts_CN %>%
    filter(chr == thisSeg$chr[1] &
             from >= thisSeg$from[1] &
             to <= thisSeg$to[1])

  ws = nrow(which_muts)

  # thisSeg
  # pio::pioStr("Chromosome", thisSeg$chr[1], suffix = '')
  # pio::pioStr(" from", thisSeg$from[1], suffix = '')
  # pio::pioStr(" to", thisSeg$to[1], suffix = '')
  # pio::pioStr(" : n =", ws, suffix = '')

  # load centromers data -- these are in absolute location format, so we update the from/ to locations
  if(offset_around_centromers > 0) {

    # Get coordinates and offset the centromers by "offset"
    data('chr_coordinates_hg19', package = 'mvmobster')

    chr_coordinates_hg19 = chr_coordinates_hg19 %>%
      mutate(
        centromerStart = centromerStart - offset_around_centromers,
        centromerEnd = centromerEnd + offset_around_centromers
      ) %>%
      filter(chr == thisSeg$chr[1])

    # Coordinates are absolute, not relative, so we need now to
    # compute the adjustment for the one that we are checking
    starts = chr_coordinates_hg19$from
    names(starts) = chr_coordinates_hg19$chr

    which_muts$from = which_muts$from + starts[1]
    which_muts$to = which_muts$to + starts[1]

    which_muts = which_muts %>%
      filter(to < chr_coordinates_hg19$centromerStart[1] |
               from > chr_coordinates_hg19$centromerEnd[1])

    ws_nc = nrow(which_muts)
    # pio::pioStr("off centromer", ws_nc, suffix = '\n')
  }

  which_muts
}

minor = function(x, seg_id, samples = x$samples)
{
  as.numeric(
    x$segments %>% filter(id == seg_id &
                            sample %in% samples &
                            variable == 'minor') %>% pull(value)
  )
}

Major = function(x, seg_id, samples = x$samples)
{
  as.numeric(
    x$segments %>% filter(id == seg_id &
                            sample %in% samples &
                            variable == 'Major') %>% pull(value)
  )
}

# VAF Adjustment for CN and purity as 1/2*CCF
vaf_adjustCN = function(v, m, M, p, mut.allele =1)
{
  CN = m+M

  0.5 * v * ((CN-2) * p + 2) / (mut.allele * p)
}
