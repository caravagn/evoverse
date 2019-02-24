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
#' @param N.min Minimum number of mutations that need to map to
#' a segment. If a segment has less than \code{N.min}, the
#' segment is rejected.
#' @param description Dataset synopsis.
#' @param offset_around_centromers When mapping mutations to
#' around segments, exclude this offset around the centromers
#' of each chromosome (reference coordinates hg19).
#'
#' @return An object that represents a dataset from class
#' \code{mbs_data}.
#'
#' @export
#'
#' @examples
#' data(example_mvmobster)
#' TODO
mobster_dataset = function(
  mutations,
  segments,
  samples,
  purity,
  description = "My multi-region MOBSTER dataset",
  N.min = 500,
  offset_around_centromers = 1e6
)
{
  pio::pioHdr("mvMOBSTER dataset")
  data = mutations

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check input format for data
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # data should be dataframe
  if(!is.data.frame(data)) stop("Data must be a dataframe")

  pioStr('Mutations', paste0('N = ', nrow(data)), suffix = '\n')

  # Columns with data
  DP.columns = paste0(samples, ".DP")
  NV.columns = paste0(samples, ".NV")
  VAF.columns = paste0(samples, ".VAF")

  Major.columns = paste0(samples, ".Major")
  minor.columns = paste0(samples, ".minor")

  all.columns = c(DP.columns, NV.columns, VAF.columns)

  if(any(!(DP.columns %in% colnames(data)))) stop("Missing DP columns in data.")
  if(any(!(NV.columns %in% colnames(data)))) stop("Missing NV columns in data.")
  if(any(!(VAF.columns %in% colnames(data)))) stop("Missing VAF columns in data.")

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Mutations must have chromosomal coordinated
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  # mutations must have locations
  if(any(!(c('chr', 'from', 'to') %in% colnames(data))))
    stop("Missing location of the annotated mutations: chr, from, to.")

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # NA values in the data are notified
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # NAs in any of the columns for values
  if(any(is.na(data[, all.columns, drop = FALSE])))
  {
    message("There are NA values in some of the entries for columns: ",
            paste(all.columns, collapse = ', '),
            '\nThese will be removed.')

    data = data[
      complete.cases(data[, all.columns, drop = FALSE]), , drop = FALSE
    ]

    pioStr('Removed NAs', paste0('N = ', nrow(data)))
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Create IDs for database
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  id_col = c('id', 'ID')

  if(!all(id_col %in% colnames(data)))
    data$id = paste0('__mut', 1:nrow(data))
  else
  {
    message("Found 'id' column in data, will use that as identifier.")

    if(id_col[2] %in% colnames(data)) data$id = data$ID
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check purity and CNA format
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(any(!(names(purity) %in% samples)))
    stop("You forgot to report purity for some samples?")

  # Copy Number must havea similar format
  if(
    !is.data.frame(segments) |
    !all(c('chr', 'from', 'to') %in% colnames(segments)) |
    !all(Major.columns %in% colnames(segments)) |
    !all(minor.columns %in% colnames(segments))
  ) {
    stop("Your segments do not look in the right format.")
  }

  tib_data = as_tibble(data)
  tib_segments = as_tibble(segments)

  # if(relative.coordinates)
  # {
  #   pio::pioStr("Switching to absolute chromosomal coordinates for", 'ref. hg19')
  #
  #   data('chr_coordinates_hg19', package = 'mvmobster')
  #
  #   starts = chr_coordinate_hg19$from
  #   names(starts) = chr_coordinate_hg19$chr
  #
  #   tib_data = tib_data %>% mutate(
  #     relative.from = from,
  #     relative.to = to,
  #     from = from + starts[chr],
  #     to = to + starts[chr]
  #   )
  #
  #   tib_segments = tib_segments %>% mutate(
  #     relative.from = from,
  #     relative.to = to,
  #     from = from + starts[chr],
  #     to = to + starts[chr]
  #   )
  # }

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Melt data
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  # Melt data
  pioStr("Melting data", "")

  tib_data = tib_data %>%
    reshape2::melt(id = 'id') %>%
    as_tibble

  cat("OK\n")

  # make everything a chr
  tib_data$variable = paste(tib_data$variable)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Create a mbst_data object
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  x = list()
  class(x) <- "mbst_data"

  x$samples = samples

  x$description = description

  # Annotations: anything but VAF/ CN values
  x$annotations = tib_data %>%
    filter(!(variable %in%
               c(all.columns, 'chr', 'from', 'to')))

  # Data: DP, NV and VAF
  x$data = tib_data %>%
    filter(variable %in% all.columns)

  x$data$value = as.numeric(x$data$value)

  sp = strsplit(x$data$variable, '\\.')
  x$data$variable = sapply(sp, function(w) w[2])
  x$data$sample = sapply(sp, function(w) w[1])

  # Locations: CN position of mutations
  x$locations = tib_data %>%
    filter(variable %in%  c('chr', 'from', 'to'))

  # Log creation
  x = logOp(x, "Initialization")

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Mapping mutations to CNAs
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #  tibble for the segments, we ID them
  tib_segments$id = paste0("_CN_seg", 1:nrow(tib_segments))

  tib_segments = tib_segments %>%
    reshape2::melt(id = c('id', 'chr', 'from', 'to')) %>%
    as_tibble

  tib_segments$variable = paste(tib_segments$variable)

  sp = strsplit(tib_segments$variable, '\\.')
  tib_segments$variable = sapply(sp, function(w) w[2])
  tib_segments$sample = sapply(sp, function(w) w[1])

  # store segments in the obj
  x$segments = tib_segments

  # store purity
  x$purity = purity

  # check that all samples have segments available
  if(!all(x$samples %in% unique(x$segments$sample)))
    stop("Some samples are missing from the list of segments, will have to stop.")

  if(any(!(unique(x$segments$sample) %in% x$samples)))
  {
    message("Your list of segments contains more sample IDs that those used -- will be dropped.")

    x$segments = x$segments %>%
      filter(sample %in% x$samples)
  }

  # we create a map for each mutation to each segment
  pio::pioTit(paste0("Mapping mutations to CN segments."))
  pio::pioStr("Offset around centromers (hg19)", offset_around_centromers, suffix = '\n', prefix = '')

  if(nrow(segments) > 1) pb = txtProgressBar(min = 0, max = nrow(segments), style = 3)
  segments_ids = unique(x$segments$id)

  x$map_mut_seg = NULL

  for(s in seq(segments_ids)) {
    if(nrow(segments) > 1) setTxtProgressBar(pb, s)

    mapped = byLoc(x, segments_ids[s], offset_around_centromers)
    if(nrow(mapped) == 0) next;

    mapped$seg_id = segments_ids[s]

    x$map_mut_seg = bind_rows(x$map_mut_seg, mapped)
  }

  if(nrow(x$map_mut_seg) == 0) {
    stop("None of the input mutations map to a segment, cannot do anything.")
  }

  # These have not been map to any segment
  unmapped = VAF(x) %>%
    filter(!id %in% x$map_mut_seg$id) %>% pull(id)
  unmapped = unique(unmapped)

  if(length(unmapped) > 0)
  {
    num = length(unmapped)
    pnum = num/N(x) * 100

    pio::pioStr(
      '\n\nOutside segments',
      paste0('N = ', num),
      suffix = paste0('(', round(pnum, 0), '%)\n')
    )

    x = delete_entries(x, unmapped)
  }

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

  return(x)
}





# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Private auxiliary functions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)

  obj
}



# Conditional mutate
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}



