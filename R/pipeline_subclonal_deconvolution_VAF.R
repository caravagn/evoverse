#' Pipeline to perform subclonal deconvolution with MOBSTER.
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
#' # We use auto_setup = 'FAST' to speed up the analysis
#' x = pipeline_subclonal_deconvolution_VAF(mutations, cna = cna, purity = purity,  N_max = 1000, auto_setup = 'FAST')
#' print(x)
pipeline_subclonal_deconvolution_VAF = function(mutations,
                                      cna = NULL,
                                      purity = NULL,
                                      min_VAF = 0.05,
                                      karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                      min_muts = 50,
                                      description = "VAF Subclonal deconvolution dataset (MOBSTER + BMix)",
                                      N_max = 15000,
                                      ...
)
{

  pio::pioHdr("Evoverse", italic(paste0('~ Subclonal deconvolution pipeline for VAF data')))
  cat('\n')

  #
  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  CNAqc_input = deconvolution_prepare_input(mutations, cna, purity, N_max, min_VAF = min_VAF)

  #
  # 2) MOBSTER analysis of karyotypes
  #
  mobster_fits = deconvolution_mobster_karyotypes_VAF(
      mutations = CNAqc_input$snvs,
      karyotypes = karyotypes,
      min_muts = min_muts,
      ...
    )

  # Extract the best fits that we are going to qc
  cat("\n")
  cli::cli_h2("MOBSTER QC")
  cat("\n")

  mobster_best_fits = lapply(
    mobster_fits,
    function(x) {
      if (x %>% is.null %>% all) return(NULL)
      x$best
      })

  # QC
  mobster_best_fits =  lapply(mobster_best_fits,
                              evoverse:::qc_deconvolution_mobster,
                              type = 'D')

  # Report a message from QC
  for(x in karyotypes)
  {
    if (mobster_best_fits[[x]] %>% is.null %>% all)
      cli::cli_alert_warning("{.field {x}}: not analysed.")
    else
    {
      if(mobster_best_fits[[x]]$QC == "PASS")
        cli::cli_alert_success("{.field {x}}: QC PASS. p = {.value {mobster_best_fits[[x]]$QC_prob}}")
      else
        cli::cli_alert_danger("{.field {red(x)}}: QC FAIL. p = {.value {mobster_best_fits[[x]]$QC_prob}}")
    }
  }

  #
  # 3) BMix -- Binomial model -- on each karyotype non-tail mutations
  #
  bmix_best_fits = deconvolution_nonread_counts_bmix_VAF(
    x = mobster_best_fits,
    karyotypes = karyotypes,
    min_muts =  min_muts
    )
  names(bmix_best_fits) = names(mobster_best_fits)

  #
  # 4) Results assembly
  #
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()

  # Fits
  results$mobster = mobster_fits
  results$bmix = bmix_best_fits

  # Input
  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input)

  # Plot and tables
  results$table$summary = deconvolution_table_summary(mobster_best_fits, bmix_best_fits)

  # Clustering assignments
  results$table$clustering_assignments = deconvolution_table_assignments(mobster_best_fits, bmix_best_fits)

  # Plot
  results$figure = deconvolution_plot_assembly(
    mobster_best_fits,
    CNAqc_input,
    bmix_best_fits,
    figure_caption = paste0("", Sys.time(), '. evoverse pipeline for subclonal deconvolution. QC: ', results$table$summary$QC_type[1]),
    figure_title = description)

  # Data id
  results$description = description

  cli::cli_process_done()

  return(results)
}

