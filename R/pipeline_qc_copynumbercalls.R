# mutations = CNAqc::example_dataset_CNAqc$snvs
# cna = CNAqc::example_dataset_CNAqc$cna
# purity = CNAqc::example_dataset_CNAqc$purity
# reference = CNAqc::example_dataset_CNAqc$reference
#
# w = evoverse:::pipeline_qc_copynumbercalls(mutations, cna, purity, smooth = F, collate = T)

#' Title
#'
#' @param mutations
#' @param cna
#' @param purity
#' @param reference
#' @param sample
#' @param matching_epsilon_peaks
#' @param output
#' @param smooth
#' @param cex
#' @param collate
#'
#' @return
#' @export
#'
#' @examples
pipeline_qc_copynumbercalls = function(
  mutations,
  cna,
  purity,
  reference = 'GRCh38',
  sample = "PAT00001",
  output = paste0(sample, '_evoverse_qcreport.pdf'),
  matching_epsilon_peaks = 0.025,
  smooth = FALSE,
  cex = .7,
  collate = FALSE)
{
  pio::pioHdr("Evoverse", italic(paste0('~ Pipeline to determine thq quality of somatic mutation and CNA data')))
  cat('\n')

  options(CNAqc_cex = cex)

  # Creating input CNAqc
  cli::cli_h1("Creating input object (smoothing {.field {smooth}}) for sample {.field {sample}}")
  cat("\n")

  x = CNAqc::init(mutations, cna, purity, ref = reference)

  print(x)

  if(smooth) x = CNAqc::smooth_segments(x)

  USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
  prevalent_karyo = CNAqc:::get_prevalent_karyotype(x, rank = 1)

  # peak detection ~ QC purity/ ploidy
  cli::cli_h1("Purity/ ploidy QC via peak detection")
  cat("\n")

  x = CNAqc::analyze_peaks(x,
                           matching_epsilon = matching_epsilon_peaks,
                           karyotypes = USE_KARYOTYPES)

  # CCF values
  cli::cli_h1("CCF estimation and QC")
  cat("\n")

  x = CNAqc::compute_CCF(x, karyotypes = USE_KARYOTYPES)

  # Fragmentation (after smoothing)
  cli::cli_h1("Detecting patterns of overfragmentation")
  cat("\n")

  x = CNAqc::detect_arm_overfragmentation(x)
  # x =  CNAqc::detect_wg_overfragmentation(x)

  # Final QC is inside CNAqc
  xqc = CNAqc:::compute_QC_table(x)

  QC_table = xqc$QC_table
  pPASS = xqc$percentage_PASS
  NA_tests = xqc$NA_tests

  # Return this information
  fit = list()
  fit$fit = x

  fit$QC$QC_table = QC_table
  fit$QC$f_PASS = pPASS
  fit$QC$NA_tests = NA_tests


  # Figure assembly
  cat("\n")
  cli::cli_h1("PDF report assembly")
  cli::cli_process_start("Creating report in file {.field {output}}")
  cat('\n')

  suppressWarnings(
    evoverse:::report_multipage_cnaqc_pipeline(
      x,
      f = output,
      cex = cex,
      sample = sample,
      collate = collate,
      score = fit$QC$f_PASS
    )
  )
  cli::cli_process_done()

  fit$plot_file = output


  return(fit)
}

get_pqc_cna_QC_table = function(x)
{
  all_karyptypes = x$peaks_analysis$fits %>% names

  if(all(is.null(x$peaks_analysis)))
  {
    peaks_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'Peaks',
      stringsAsFactors = FALSE
    )
  }
  else
  {
    peaks_QC = x$peaks_analysis$matches %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
      dplyr::distinct(karyotype, QC, .keep_all = T) %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(
        value = 1,
        lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        label = karyotype,
        type = 'Peaks')
  }

  if(all(is.null(x$CCF_estimates)))
  {
    CCF_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'CCF',
      stringsAsFactors = FALSE
    )
  }
  else
  {
    CCF_QC = Reduce(dplyr::bind_rows, lapply(x$CCF_estimates, function(x) x$QC_table)) %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(
        value = 1,
        lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        label = karyotype,
        type = 'CCF')
  }

  QC_table = dplyr::bind_rows(peaks_QC, CCF_QC)
  QC_table$karyotype = factor(QC_table$karyotype, all_karyptypes)
  QC_table$type = factor(QC_table$type, levels = c('Peaks', 'CCF'))

  QC_table %>%
    dplyr::mutate(
      QC = ifelse(!is.na(QC) & QC == "NA", NA, QC)
    )
}

