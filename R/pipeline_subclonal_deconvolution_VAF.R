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
#' @param min_muts Skip analysing karyotypes with less than \code{min_muts} mutations.
#' @param ... Parameters passed \code{mobster_fit} in package \code{mobster}.
#' @param sample
#' @param output
#' @param min_VAF
#' @param karyotypes
#' @param N_max
#' @param matching_epsilon_peaks
#' @param cex
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
                                                sample = "XXX_YYY_ZZZ",
                                                output = paste0(sample, '_evoverse_vaf_deconvolution.pdf'),
                                                min_VAF = 0.05,
                                                karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                                min_muts = 50,
                                                N_max = 15000,
                                                matching_epsilon_peaks = 0.025,
                                                cex = .7,
                                                collate = FALSE,
                                                ...)
{
  pio::pioHdr("Evoverse", italic(paste0(
    '~ Subclonal deconvolution pipeline for VAF data'
  )))
  cat('\n')

  # Set global cex value for all next plots
  options(CNAqc_cex = cex)
  options(mobster_cex = cex)
  options(bmix_cex = cex)


  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  CNAqc_input = evoverse:::deconvolution_prepare_input(mutations, cna, purity, N_max, min_VAF = min_VAF)

  # Determine karyotypes with sufficiently good quality data
  CNAqc_input = CNAqc::analyze_peaks(CNAqc_input, matching_epsilon = matching_epsilon_peaks)

  # 2) MOBSTER analysis of karyotypes ~ all analysed
  mobster_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
    mutations = CNAqc_input$snvs,
    karyotypes = karyotypes,
    min_muts = min_muts,
    QC_type = "D",
    ...
  )

  # 3) BMix -- Binomial model -- on each karyotype, non-tail mutations
  bmix_best_fits = evoverse:::deconvolution_nonread_counts_bmix_VAF(x = mobster_fits,
                                                                    karyotypes = karyotypes,
                                                                    min_muts =  min_muts)
  names(bmix_best_fits) = names(mobster_fits)

  # 4) Results assembly and final QC
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()
  results$type = "VAF pipeline"
  class(results) = "evopipe_deconv"

  # Fits
  results$mobster = mobster_fits
  results$bmix = bmix_best_fits

  # Input
  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input
  )

  # Plot and tables
  results$table$summary = evoverse:::deconvolution_table_summary(mobster_fits, bmix_best_fits)

  # Clustering assignments
  results$table$clustering_assignments = evoverse:::deconvolution_table_assignments(mobster_fits, bmix_best_fits)

  # Data id
  results$description = sample
  results$log = paste0(
    "",
    Sys.time(),
    '. evoverse pipeline for subclonal deconvolution from VAF data. QC: ',
    results$table$summary$QC_type[1]
  )

  ##### ##### ##### #####
  ##### Final QC
  ##### ##### ##### #####

  # We want to take a decision based only on the karyotypes that pass peak detection
  QC_peaks_karyotypes = CNAqc_input$peaks_analysis$matches %>%
    dplyr::distinct(karyotype, QC) %>%
    dplyr::rename(QC_peaks = QC)

  # We want to include the karyotype size in the weighted decision
  QC_karyotypes_size = data.frame(
    karyotype = CNAqc_input$n_karyotype %>% unlist %>% names,
    N = CNAqc_input$n_karyotype,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(karyotype %in% karyotypes) %>%
    dplyr::rename(N = N.Freq) %>%
    dplyr::select(-N.Var1)

  # We want to decide only if mobster passes the fit
  QC_mobster_karyotypes = results$table$summary %>%
    dplyr::mutate(
      architecture =
        case_when(
          karyotype %in% c("1:0", "1:1") & BMix_K_B  > 1 ~ "Polyclonal",
          karyotype %in% c("1:0", "1:1") &
            BMix_K_B  == 1 ~ "Monoclonal",
          karyotype %in% c("2:0", "2:1", "2:2") &
            BMix_K_B  > 2 ~ "Polyclonal",
          karyotype %in% c("2:0", "2:1", "2:2") &
            BMix_K_B  <= 2 ~ "Monoclonal"
        )
    ) %>%
    dplyr::select(karyotype,
                  QC,
                  QC_type,
                  QC_prob,
                  architecture,
                  K_beta,
                  BMix_K_B) %>%
    dplyr::rename(K_mobster = K_beta, K_BMix = BMix_K_B)

  QC_table = data.frame(karyotype = karyotypes, stringsAsFactors = F) %>%
    dplyr::left_join(QC_peaks_karyotypes, by = 'karyotype') %>%
    dplyr::left_join(QC_mobster_karyotypes, by = 'karyotype') %>%
    dplyr::left_join(QC_karyotypes_size, by = 'karyotype')

  # Marginalise the architecture,
  final_architecture = QC_table %>%
    dplyr::filter(QC_peaks == "PASS") %>%
    dplyr::group_by(architecture) %>%
    dplyr::summarise(score = sum(N)) %>%
    dplyr::arrange(desc(score))

  # Apply pplicy
  decided_final_architecture = NA
  QC_final = 'FAIL'

  if (nrow(final_architecture) > 0)
  {
    # Largest one
    decided_final_architecture = final_architecture$architecture[1]

    # Take the largest
    QC_final = QC_table %>%
      dplyr::filter(architecture == decided_final_architecture) %>%
      dplyr::arrange(desc(N))

    QC_final = (QC_final$QC[1] == QC_final$QC_peaks) &
      (QC_final$QC[1] == "PASS")
    QC_final = ifelse(QC_final, "PASS", "FAIL")
  }

  # Store all
  results$table$QC_table = QC_table %>% as_tibble()
  results$table$QC = QC_final
  results$table$architecture = decided_final_architecture

  # Reports
  cat("\n")
  cli::cli_h2("Summary")

  pio::pioDisp(QC_table)

  cat('\n')
  cli::cli_h3("Tumour architecture: {.field {decided_final_architecture}}.")

  evoverse:::report_multipage_deconv_pipeline(
    x = results,
    f = output,
    sample = sample,
    collate = collate
  )

  cli::cli_process_done()


  return(results)
}
