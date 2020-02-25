# setwd("~/Documents/Davros/Golden")
# x = 'Golden_CS18_T'

#' Parese Sequenza Copy Number data.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
evoparse_Sequenza_CNAs = function(x)
{
  cli::cli_process_start("Evoverse parser for Sequenza CNA calls for sample {.field {x}}")
  cat('\n')

  # Files required
  segments_file = paste0(x, "_segments.txt")
  mutations_file = paste0(x, "_mutations.txt")
  purity_file = paste0(x, "_confints_CP.txt")
  alternatives_file = paste0(x, "_alternative_solutions.txt")

  if(!file.exists(segments_file))
    stop(paste0("Missing required segments file: ", segments_file))

  if(!file.exists(mutations_file))
    stop(paste0("Missing required mutations file: ", mutations_file))

  if(!file.exists(purity_file))
    stop(paste0("Missing required solutions file: ", purity_file))

  if(!file.exists(alternatives_file))
    stop(paste0("Missing alternative solutions file: ", alternatives_file))

  # Segments data -- with CNAq format conversion which is used in evoverse
  segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
    dplyr::rename(
      chr = chromosome,
      from = start.pos,
      to = end.pos,
      Major = A,
      minor = B
      ) %>%
    dplyr::select(chr, from, to, Major, minor, dplyr::everything())

  # Mutations data -- with CNAq format conversion which is used in evoverse
  mutations = readr::read_tsv(mutations_file, col_types = readr::cols()) %>%
    dplyr::rename(
      chr = chromosome,
      from = position,
    ) %>%
    dplyr::rowwise() %>%
    tidyr::separate(mutation, sep = '>', into = c('ref', 'alt'), remove = FALSE) %>%
    dplyr::mutate(to = from + nchar(alt)) %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything())

  # Solution data
  solutions = readr::read_tsv(purity_file, col_types = readr::cols())

  purity = solutions$cellularity[2]
  ploidy = solutions$ploidy.estimate[2]

  # Alternative solutions data
  alternative_solutions = readr::read_tsv(alternatives_file, col_types = readr::cols())

  fits = list()
  fits$input = x
  fits$segments = segments
  fits$mutations = mutations
  fits$purity = purity
  fits$ploidy = ploidy
  fits$alternative_solutions = alternative_solutions

  cli::cli_process_done()

  return(fits)
}

#' Pares Platypus mutation data.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
evoparse_platypus_mutations = function(x)
{
  cli::cli_process_start("Evoverse parser for Platypus mutation calls from VCF file {.field {x}}")
  cat('\n')

  if(!file.exists(x))
    stop(paste0("Missing required VCF file: ", x))

  # Parsing with vcfR -- maybe Variant Annontation works better..
  xfile = vcfR::read.vcfR(x)

  # INFO field does not come out very nice -- so I omit it
  info_field = vcfR::extract_info_tidy(xfile)

  # Genotype(s)
  gt_field = vcfR::extract_gt_tidy(xfile) %>%
    dplyr::mutate(
      DP = as.numeric(gt_NR),
      NV = as.numeric(gt_NV),
      VAF = NV/DP
    )

  na_muts = gt_field %>% dplyr::filter(is.na(VAF))

  if(nrow(na_muts) > 0) {
    cat("There are n =", nrow(na_muts), "mutations with VAF that is NA; they will NOT be removed.\n")
  }

  # Fixed fields -- this should retain everything relevant
  fix_field = xfile@fix %>%
    as_tibble() %>%
    dplyr::rename(
      chr = CHROM,
      from = POS,
      ref = REF,
      alt = ALT
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(alt)
      ) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything())

  if(nrow(fix_field) != nrow(gt_field))
    stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

  fits = list()
  fits$input = x
  fits$mutations = dplyr::bind_cols(fix_field, gt_field) %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())

  cli::cli_process_done()

  return(fits)
}

# Sequenza_calls = evoparse_Sequenza_CNA(x)


