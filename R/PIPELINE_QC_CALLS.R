#
# mutations = CNAqc::example_dataset_CNAqc$snvs
# cna = CNAqc::example_dataset_CNAqc$cna
# purity = CNAqc::example_dataset_CNAqc$purity
# reference = CNAqc::example_dataset_CNAqc$reference
#
# w = evoverse:::pipeline_qc_copynumbercalls(mutations, cna, purity, smooth = T)
# evoverse:::pdf_report(w, file = 'a.pdf')
# plot(w, cex = 1)

#' Pipeline for QC of somatic calls (copy number and mutations).
#'
#' @description
#'
#' This is ...
#'
#' All input data must be in the format supported by package \code{CNAqc}.
#'
#' @param mutations Somatic mutation calls, in \code{CNAqc} input format.
#' @param cna Absolute Copy Number segments, in \code{CNAqc} input format.
#' @param purity Tumour purity, in \code{CNAqc} input format.
#' @param reference A reference supported by package \code{CNAqc}.
#' @param description A sample description; will appear in the plots.
#' @param purity_error Error on purity accepted \code{CNAqc}.
#' @param VAF_tolerance Accepted tolerance on VAF picks calling \code{CNAqc}.
#' @param smooth If \code{TRUE}, smooths copy number data with \code{CNAqc}.
#' @param ccf_method
#' @param peak_method
#' @param min_CCF
#' @param only_SNVs
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
  description = "My sample",
  smooth = TRUE,
  purity_error = 0.05,
  VAF_tolerance = 0.015,
  ccf_method = 'ENTROPY',
  peak_method = 'closest',
  only_SNVs = TRUE
 )
{
  pio::pioHdr("Evoverse", crayon::italic(paste0('~ Pipeline to QC somatic mutations and CNA segments')))
  cat('\n')

  # 1. Creating input CNAqc
  cli::cli_h1("Creating input object (smoothing {.field {smooth}}) for sample {.field {description}}")
  cat("\n")

  x = evoverse:::deconvolution_prepare_input(mutations, cna, purity, reference, only_SNVs = only_SNVs)
  print(x)

  if(smooth) x = CNAqc::smooth_segments(x)

  USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
  prevalent_karyo = CNAqc:::get_prevalent_karyotype(x, rank = 1)

  # 2. peak detection ~ QC purity/ ploidy
  cli::cli_h1("Purity/ ploidy QC via peak detection")
  cat("\n")

  x = CNAqc::analyze_peaks(
    x,
    purity_error = purity_error,
    VAF_tolerance = VAF_tolerance, 
    karyotypes = USE_KARYOTYPES,
    matching_strategy = peak_method
    )

  # 3. CCF values
  cli::cli_h1("CCF estimation and QC")
  cat("\n")

  x = CNAqc::compute_CCF(x, karyotypes = USE_KARYOTYPES, method = ccf_method)

  # 4. Fragmentation (after smoothing)
  cli::cli_h1("Detecting patterns of overfragmentation")
  cat("\n")

  x = CNAqc::detect_arm_overfragmentation(x)
  # x =  CNAqc::detect_wg_overfragmentation(x)
  
  # 5. Advanced phasing
  cli::cli_h1("Calculating multiplicieties for complex karyotypes")
  phas = CNAqc:::advanced_phasing(x)
  x$advanced_phasing = phas$phasing

  # 6. Final QC is inside CNAqc
  xqc = CNAqc:::compute_QC_table(x)

  # Assemble a special object to return this information
  fit = list()
  class(fit) = "evopipe_qc"

  fit$type = "QC pipeline"
  fit$description = description

  fit$cnaqc = x

  fit$QC$QC_table = xqc$QC_table
  fit$QC$f_PASS = xqc$percentage_PASS
  fit$QC$NA_tests = xqc$NA_tests

  fit$log = paste0("", Sys.time(), '. evoverse pipeline for data QC.')

  return(fit)
}

#' S3 print for the pipeline results
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.evopipe_qc = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_qc'))

  # Print pipeline objects
  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ Evoverse ] {.value {x$description}}")),
      '{.value {x$type}}'
    )
  )
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$cnaqc)

  # PASS/FAIL
  pass = x$QC$QC_table %>% dplyr::filter(QC == "PASS")
  fail = x$QC$QC_table %>% dplyr::filter(QC == "FAIL")

  if(nrow(pass) > 0) {
    cat("\n")
    cli::cli_rule(crayon::bgGreen(" QC PASS "), right = paste0("PASS rate (%): ", x$QC$f_PASS))
    print(pass)
  }

  if(nrow(fail) > 0) {
    cat("\n")
    cli::cli_rule(crayon::bgRed(" QC PASS "), right = paste0("FAIL rate (%): ", 100 - x$QC$f_PASS))
    print(fail)
  }
}

# S3 plot for the pipeline results
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.evopipe_qc = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_qc'))

  # Segment plots
  segments_plot = cowplot::plot_grid(
    suppressWarnings(suppressMessages(CNAqc::plot_gw_depth(x$cnaqc))),
    suppressWarnings(suppressMessages(CNAqc::plot_segments(x$cnaqc, highlight = c('1:0', '1:1', '2:0', '2:1', '2:2')))),
    align = 'v',
    nrow = 2,
    rel_heights = c(.15, .8)
  )

  # Data histograms
  hist_plot = ggpubr::ggarrange(
    suppressWarnings(suppressMessages(CNAqc::plot_data_histogram(x$cnaqc, which = 'CCF'))),
    suppressWarnings(suppressMessages(CNAqc::plot_data_histogram(x$cnaqc, which = 'VAF'))),
    suppressWarnings(suppressMessages(CNAqc::plot_data_histogram(x$cnaqc, which = 'DP'))),
    suppressWarnings(suppressMessages(CNAqc::plot_data_histogram(x$cnaqc, which = 'NV'))),
    nrow = 1,
    ncol = 4,
    common.legend = TRUE,
    legend = 'bottom'
  )

  peak_plot =  suppressWarnings(suppressMessages(CNAqc::plot_peaks_analysis(x$cnaqc)))
  CCF_plot =  suppressWarnings(suppressMessages(CNAqc::plot_CCF(x$cnaqc, strip = TRUE)))

  # One page
  figure = ggpubr::ggarrange(
    segments_plot,
    hist_plot,
    peak_plot,
    CCF_plot,
    ncol = 1,
    nrow = 4,
    heights = c(1.1, .9, .9, .9)
  )

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(
    figure,
    bottom = ggpubr::text_grob(
      bquote(.(x$log)),
      hjust = 0,
      x = 0,
      size = 10
    )
  )

  return(figure)
}


