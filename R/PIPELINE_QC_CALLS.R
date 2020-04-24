#
# mutations = CNAqc::example_dataset_CNAqc$snvs
# cna = CNAqc::example_dataset_CNAqc$cna
# purity = CNAqc::example_dataset_CNAqc$purity
# reference = CNAqc::example_dataset_CNAqc$reference
#
# w = evoverse:::pipeline_qc_copynumbercalls(mutations, cna, purity, smooth = T)
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
#' @param sample A sample identifier; will appear in the plots.
#' @param matching_epsilon_peaks Epsilon to match peaks with \code{CNAqc}.
#' @param smooth If \code{TRUE}, smooths copy number data with \code{CNAqc}.
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
  sample = "MyPAT00XX1",
  smooth = TRUE,
  matching_epsilon_peaks = 0.025,
  CCF_computation = 'ENTROPY'
 )
{
  pio::pioHdr("Evoverse", crayon::italic(paste0('~ Pipeline to QC somatic mutations and CNA segments')))
  cat('\n')

  # 1. Creating input CNAqc
  cli::cli_h1("Creating input object (smoothing {.field {smooth}}) for sample {.field {sample}}")
  cat("\n")

  x = CNAqc::init(mutations, cna, purity, ref = reference)

  print(x)

  if(smooth) x = CNAqc::smooth_segments(x)

  USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
  prevalent_karyo = CNAqc:::get_prevalent_karyotype(x, rank = 1)

  # 2. peak detection ~ QC purity/ ploidy
  cli::cli_h1("Purity/ ploidy QC via peak detection")
  cat("\n")

  x = CNAqc::analyze_peaks(
    x,
    matching_epsilon = matching_epsilon_peaks,
    karyotypes = USE_KARYOTYPES
    )

  # 3. CCF values
  cli::cli_h1("CCF estimation and QC")
  cat("\n")

  x = CNAqc::compute_CCF(x, karyotypes = USE_KARYOTYPES, method = CCF_computation)

  # 4. Fragmentation (after smoothing)
  cli::cli_h1("Detecting patterns of overfragmentation")
  cat("\n")

  x = CNAqc::detect_arm_overfragmentation(x)
  # x =  CNAqc::detect_wg_overfragmentation(x)

  # 5. Final QC is inside CNAqc
  xqc = CNAqc:::compute_QC_table(x)

  # Assemble a special object to return this information
  fit = list()
  class(fit) = "evopipe_qc"

  fit$type = "QC pipeline"
  fit$description = sample

  fit$fit = x

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
  y = x$fit
  y$peaks_analysis = y$CCF_estimates = y$arm_fragmentation = NULL
  CNAqc:::print.cnaqc(y)

  # PASS/FAIL
  pass = x$QC$QC_table %>% dplyr::filter(QC == "PASS") %>% dplyr::select(-lab.ypos, -label, -value)
  fail = x$QC$QC_table %>% dplyr::filter(QC == "FAIL") %>% dplyr::select(-lab.ypos, -label, -value)

  if(nrow(pass) > 0) {
    cli::cli_rule(crayon::bgGreen(" QC PASS "), right = paste0("PASS rate (%): ", x$QC$f_PASS))
    print(pass)
  }

  if(nrow(fail) > 0) {
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
  get_param_ell = function(x, ...) {
    p = list(...)
    stopifnot(x %in% names(p))
    p[x]
  }

  stop("Does Not work no save to file")

  stopifnot(inherits(x, 'evopipe_qc'))

  # Figure assembly
    suppressWarnings(
      evoverse:::report_multipage_cnaqc_pipeline(
        x,
        f = file,
        cex = get_param_ell(x = 'cex', ...),
        sample = x$description,
        collate = TRUE,
        score = x$QC$f_PASS
      )
    )

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wrapper function for evopipe_qc
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cli::cli_process_done()

}


