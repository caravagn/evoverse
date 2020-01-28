#
# cna = CNAqc::example_dataset_CNAqc$cna
# mutations = CNAqc::example_dataset_CNAqc$snvs
# purity = CNAqc::example_dataset_CNAqc$purity

# x = chromosome_timing_pipeline(mutations, cna = cna, purity = purity, auto_setup = 'FAST')

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
#' x = pipeline_chromosome_timing(mutations, cna = cna, purity = purity, auto_setup = 'FAST')
#' print(x)
pipeline_chromosome_timing = function(mutations,
                                      cna = NULL,
                                      purity = NULL,
                                      timeable = c('2:0', '2:1', '2:2'),
                                      min_muts = 50,
                                      ...
                                      )
{
  pio::pioHdr("Evoverse", italic('Copy Number timing pipeline'))
  cat('\n')

  # x must be a dataset for MOBSTER
  if(!is.data.frame(mutations) & !is.matrix(mutations))
    stop("mutations: Must use a dataset for MOBSTER!")

  if(!is.null(cna) & !is.data.frame(cna) & !is.matrix(cna))
    stop("cna: Must use a dataset for MOBSTER!")

  if(!is.null(purity) & (purity > 1 | purity <= 0))
    stop("Purity must be in 0/1.")

  # Apply CNA mapping and retain only mappable mutations
  if(!is.null(cna))
  {
    cli::cli_process_start("Using CNA data to subset mutations.")
    cat("\n")

    cna_obj = CNAqc::init(mutations, cna, purity)
    mutations = cna_obj$snvs

    cli::cli_process_done()
  }

  # Deconvolution function
  mfits = deconvolution_mobster_karyotypes(
    mutations = mutations,
    karyotypes = timeable,
    min_muts = min_muts,
    ...
  )

  # QC with the trained classifier evoverse::qc_timing_model (type = 'T')
  cat("\n")
  cli::cli_rule("QC MOBSTER fits results")

  qc_clocks =  lapply(mfits, function(d) qc_deconvolution_mobster(d$best, type = 'T'))

  qc_clocks_table = lapply(qc_clocks,
                           function(x)
                             bind_cols(
                               mobster::to_string(x),
                               data.frame(QC = x$QC, QC_prob = x$QC_prob, QC_type = x$QC_type)
                               )
                           )
  qc_clocks_table = Reduce(bind_rows, qc_clocks_table)

  for(l in seq(qc_clocks_table$QC))
  {
    if(qc_clocks_table$QC[l] == "PASS")
      cli::cli_alert_success("Karyotype {.field {qc_clocks_table$karyotype[l]}} QC PASS. p = {.value {qc_clocks_table$QC_prob[l]}}")
    else
      cli::cli_alert_danger("Karyotype {.field {red(qc_clocks_table$karyotype[l])}} QC FAIL. p = {.value {qc_clocks_table$QC_prob[l]}}")
  }

  # Produce plots
  cli::cli_process_start("Preparing plots and tables for MOBSTER fits.")

  mfits_plot = lapply(qc_clocks, function(y) qc_mobster_plot(y))

  cna_plot = ggplot() + geom_blank()
  timeable = c("WG", timeable)

  if(!is.null(cna)) cna_plot = CNAqc::plot_icon_CNA(cna_obj)

  mfits_plot =  append(cna_plot, list(mfits_plot))
  mfits_plot = ggpubr::ggarrange(plotlist = mfits_plot,
                                 nrow = 1, ncol = length(mfits_plot), labels = timeable)


  # Table of assignments
  atab =  Reduce(bind_rows,
                 lapply(mfits,
                        function(x) {
                          if (all(is.null(x)))
                            return(NULL)

                          mobster::Clusters(x$best)
                        }))

  cli::cli_process_done()

  return(
    list(
      input = list(mutations = mutations, cna = cna, purity = purity),
      mobster = list(fits = mfits, plots = mfits_plot),
      assignments = atab,
      qc_clocks = qc_clocks
      )
    )
}

