#' Pipeline to time aneuploidy with MOBSTER.
#'
#' @description
#'
#' For a list of input mutations, CNA and purity data from a tumour sample,
#' MOBSTER fits are computed to determine proportions of mutations in each
#' aneuploidy state (copy state). Assignments and plots are returned, and
#' an automatic QC pipeline (supervised logistic classifier) is used to
#' flag which fits pass quality check. The trained set for the classifier
#' consists in n = 377 colorectal whole genome samples.
#'
#' @param mutations Input mutations in the format for package \code{CNAqc}.
#' @param cna Absolute Copy Number Alterations in the format for package \code{CNAqc}.
#' @param purity Tumour purity, in [0, 1].
#' @param timeable Karyotypes that can be timed, expressed as "Major:minor" notation. Default
#' are the following states \code{timeable = c('2:0', '2:1', '2:2')}.
#' @param min_muts Skip analysing karyotypes with less than \code{min_muts} mutations.
#' @param ... Parameters passed \code{mobster_fit} in package \code{mobster}.
#'
#' @return
#' @export
#'
#' @examples
#' # We use data released with the CNAqc package
#'
#' cna = CNAqc::example_dataset_CNAqc$cna
#' mutations = CNAqc::example_dataset_CNAqc$snvs
#' purity = CNAqc::example_dataset_CNAqc$purity
#'
#' x = pipeline_chromosome_timing(mutations, cna = cna, purity = purity, auto_setup = 'FAST', N_max = 500)
#' print(x)
pipeline_chromosome_timing = function(mutations,
                                      cna,
                                      purity,
                                      timeable = c('2:0', '2:1', '2:2'),
                                      min_muts = 50,
                                      min_VAF = 0.05,
                                      description = "Chromosomal timing dataset (MOBSTER)",
                                      N_max = 15000,
                                      ...
                                      )
{
  pio::pioHdr("Evoverse", italic('Chromosomal timing pipeline'))
  cat('\n')

  #
  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  CNAqc_input = deconvolution_prepare_input(mutations, cna, purity, N_max, min_VAF = min_VAF)

  #
  # 2) MOBSTER analysis of karyotypes, this returns only the best fit
  #
  mobster_fits = deconvolution_mobster_karyotypes_VAF(
    mutations = CNAqc_input$snvs,
    karyotypes = timeable,
    min_muts = min_muts,
    QC_type = "T",
    ...
  )

  #
  # 3) Results assembly
  #
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()
  results$type = "Timing pipeline"
  class(results) = "evopipe_deconv"

  # Fits
  results$mobster = mobster_fits

  # Input
  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input)

  # Plot and tables
  results$table$summary = deconvolution_table_summary(mobster_fits, bmix_fits = NULL)

  # Clustering assignments
  results$table$clustering_assignments = deconvolution_table_assignments(mobster_fits, bmix_fits = NULL)

  # Data id
  results$description = description
  results$log = paste0("", Sys.time(), '. evoverse pipeline for chromosome timing. QC: ', results$table$summary$QC_type[1])

  cli::cli_process_done()


  return(results)
}

