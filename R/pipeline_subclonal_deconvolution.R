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
pipeline_subclonal_deconvolution = function(mutations,
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

  # Timeable ones, processed one by one with MOBSTER
  # VAF required in 0/1
  mfits = lapply(
    timeable,
    function(k)
    {
      cat("\n")
      cli::cli_process_start("MOBSTER clustering of mutations with karyotype {.field {k}}")
      cat("\n")

      kmuts = mutations %>%
        filter(karyotype == k, VAF > 0, VAF < 1)

      if(nrow(kmuts) < min_muts) {
        cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
        cli::cli_process_failed()

        return(NULL)
      }

      mfit = mobster::mobster_fit(kmuts, ...)

      cat("\n")
      cli::cli_process_done()

      return(mfit)
    }
  )
  names(mfits) = timeable

  # QC with the trained classifier evoverse::qc_timing_model
  cat("\n")
  cli::cli_rule("QC MOBSTER fits results")

  qc_timing_model = evoverse::qc_timing_model

  qc_clocks = lapply(
    lapply(timeable, function(x) mfits[[x]]$best),
    qc_mobster_timing,
    input = input
  )

  qc_clocks = Reduce(bind_rows, qc_clocks)

  qc_clocks$karyotype = timeable

  # Results: class and probability
  class_of = pio:::nmfy(qc_clocks$karyotype, qc_clocks$QC)

  p_of = pio:::nmfy(qc_clocks$karyotype, qc_clocks$QC_prob)
  p_of = round(p_of, 3)

  # Produce plots
  cli::cli_process_start("Preparing plots and tables for MOBSTER fits.")

  mfits_plot = lapply(names(mfits),
                      function(y) {
                        if(all(is.null(mfits[[y]]))) return(ggplot() + geom_blank())

                        qc = ifelse(class_of[y] == "FAIL",
                                    "indianred3",
                                    'forestgreen')

                        mobster::plot.dbpmm(mfits[[y]]$best) +
                          labs(title = bquote("AUTO QC "~ .(class_of[y]) ~ p["PASS"] ~'='~ .(p_of[y]))) +
                          theme(
                            title = element_text(color = qc),
                            panel.border = element_rect(colour = qc, fill = NA, size = 5)
                          )

                      })

  if(!is.null(cna)) {
    mfits_plot =  append(CNAqc::plot_icon_CNA(cna_obj), list(mfits_plot))
    timeable = c("WG", timeable)
  }

  mfits_plot = ggpubr::ggarrange(plotlist = mfits_plot,
                                 nrow = 1, ncol = length(mfits_plot), labels = timeable)

  # Table of assignments
  atab =  lapply(mfits,
                 function(x) {
                   if(all(is.null(x))) return(NULL)

                   mobster::Clusters(x$best)
                 })

  atab = Reduce(bind_rows, atab)

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

qc_mobster_timing = function(x, input, model)
{
  m = x

  # Fail
  if(m$Kbeta == 1) return(0)

  scores = data.frame(
    tailTRUE = m$fit.tail,
    reduced.entropy = m$scores$reduced.entropy,
    entropy = m$scores$entropy
  )

  input = bind_cols(mobster:::to_string(m), scores)  %>%
    mutate(ratioM = abs((Mean_C1 / Mean_C2) - 2),
           ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
           ratioN = (max(N_C1, N_C2) +1) / (min(N_C2, N_C1) + 1),
           minpi = min(pi_C1, pi_C2))

  fulldatasetprob <- model %>% predict(input, type = "response")
  predicted.classes.fulldataset <- ifelse(fulldatasetprob > 0.5, "PASS", "FAIL")

  if(predicted.classes.fulldataset == "PASS")
    cli::cli_alert_success("Karyotype {.field {k}} QC PASS.")
  else
    cli::cli_alert_danger("Karyotype {.field {k}} QC FAIL.")


  input$QC <- predicted.classes.fulldataset
  input$QC_prob <- fulldatasetprob

  return(input)
}
