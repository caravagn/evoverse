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
#' x = pipeline_subclonal_deconvolution(mutations, cna = cna, purity = purity,  N_max = 1000, auto_setup = 'FAST')
#' print(x)
pipeline_subclonal_deconvolution_CCF = function(mutations,
                                            cna = NULL,
                                            purity = NULL,
                                            min_VAF = 0.05,
                                            CCF_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                            min_muts = 50,
                                            description = "Subclonal deconvolution dataset (MOBSTER + BMix)",
                                            N_max = 15000,
                                            ...
)
{

  pio::pioHdr("Evoverse", italic(paste0('~ Subclonal deconvolution pipeline for CCF data')))
  cat('\n')

  #
  # 1) Check input data, subset by CNA data if available
  #
  cat("\n")
  cli::cli_process_start("Loading input data")
  cat("\n")

  CNAqc_input = deconvolution_prepare_input(mutations, cna, purity, min_VAF = min_VAF)

  # Downsample data if too many mutations
  if(CNAqc_input$n_snvs > N_max) {

    cli::boxx(paste0("n = ", CNAqc_input$n_snvs,  " mutations. n > ", N_max, " , downsampling input.")) %>% cat
    cat('\n')

    CNAqc_input = CNAqc_input %>% CNAqc::subsample(N = N_max, keep_drivers = TRUE)
  }

  print(CNAqc_input)
  cli::cli_process_done()

  #
  # 2) MOBSTER analysis of CCF special function that computes CCF and adjust it for the VAF etc.
  #
  CCF_fit = deconvolution_mobster_CCF(CNAqc_input,
                                      CCF_karyotypes = CCF_karyotypes,
                                      min_muts = min_muts,
                                      min_VAF = min_VAF,
                                      QC_type = "D",
                                      ...)

  # Update the object that now has also CCF data
  CNAqc_input = CCF_fit$cna_obj
  mobster_fits = list(`CCF` = CCF_fit$fits)

  #
  # 3) BMix -- Binomial model -- with special adjustment for CCF
  #
  bmix_best_fits = deconvolution_nonread_counts_bmix_CCF(
    x = mobster_fits,
    min_muts =  min_muts
  )
  bmix_best_fits = list(`CCF` = bmix_best_fits)
  names(bmix_best_fits) = "CCF"

  #
  # 4) Results assembly
  #
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()
  results$type = "CCF pipeline"
  class(results) = "evopipe_deconv"

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
  results$table$summary = deconvolution_table_summary(mobster_fits, bmix_best_fits)

  # Clustering assignments
  results$table$clustering_assignments = deconvolution_table_assignments(mobster_fits, bmix_best_fits)

  # Plot
  # results$figure = deconvolution_plot_assembly(
  #   mobster_best_fits,
  #   CNAqc_input,
  #   bmix_best_fits,
  #   figure_caption = paste0("", Sys.time(), '. evoverse pipeline for subclonal deconvolution. QC: ', results$table$summary$QC_type[1]),
  #   figure_title = description)

  # Data id
  results$description = description
  results$log = paste0("", Sys.time(), '. evoverse pipeline for subclonal deconvolution from CCF data. QC: ', results$table$summary$QC_type[1])

  cli::cli_process_done()

  return(results)
}

