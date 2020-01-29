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
                                      cna = NULL,
                                      purity = NULL,
                                      timeable = c('2:0', '2:1', '2:2'),
                                      min_muts = 50,
                                      description = "Chromosomal timing dataset (MOBSTER)",
                                      N_max = 15000,
                                      ...
                                      )
{
  pio::pioHdr("Evoverse", italic('Chromosomal timing pipeline'))
  cat('\n')

  #
  # 1) Check input data, subset by CNA data if available
  #
  prepared_input = deconvolution_prepare_input(mutations, cna, purity, N_max)
  mutations = prepared_input$mutations
  cna = prepared_input$cna
  cna_obj = prepared_input$cna_obj
  purity = prepared_input$purity

  #
  # 2) MOBSTER analysis of karyotypes
  #
  mobster_fits = evoverse:::deconvolution_mobster_karyotypes(
    mutations = mutations,
    karyotypes = timeable,
    min_muts = min_muts,
    ...
  )

  # Extract the best fits that we are going to qc
  mobster_best_fits = lapply(mobster_fits,
                             function(x){
                               if (x %>% is.null %>% all)
                                 return(NULL)
                               x$best
                             })

  cat("\n")
  cli::cli_h2("MOBSTER QC")
  cat("\n")

  # QC
  mobster_best_fits =  lapply(mobster_best_fits, evoverse:::qc_deconvolution_mobster, type = 'T')

  # Report a message from QC
  for(x in timeable)
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
  # 4) Results assembly
  #
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()

  # Fits
  results$mobster = mobster_fits

  # Input
  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = cna_obj)

  # Plot and tables
  results$table$summary = deconvolution_table_summary(mobster_best_fits, bmix_fits = NULL)

  # Clustering assignments
  results$table$clustering_assignments = deconvolution_table_assignments(mobster_best_fits, bmix_fits = NULL)

  results$figure = deconvolution_plot_assembly(
    mobster_fits = mobster_best_fits,
    cna_obj,
    bmix_fits = NULL,
    figure_caption = paste0("", Sys.time(), '. evoverse pipeline for chromosome timing. QC: ', results$table$summary$QC_type[1]),
    figure_title = description)

  # Data id
  results$description = description

  cli::cli_process_done()


  return(results)
}

