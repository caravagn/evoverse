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
  pio::pioHdr("Evoverse", crayon::italic('Raw VAF by karyotype deconvolution pipeline'))
  cat('\n')

  # # Set global cex value for all next plots
  # options(CNAqc_cex = cex)
  # options(mobster_cex = cex)
  # options(bmix_cex = cex)


  #
  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  cat("\n")
  cli::cli_process_start("Loading input data")
  cat("\n")

  CNAqc_input = evoverse:::deconvolution_prepare_input(mutations,
                                                       cna,
                                                       purity,
                                                       reference = reference,
                                                       min_VAF = min_VAF)

  print(CNAqc_input)
  cli::cli_process_done()

  #
  # 2) MOBSTER + BMix analysis of karyotypes, this returns only the best fit
  #
  all_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
    x = CNAqc_input,
    karyotypes = karyotypes,    # Required karyotypes
    BMix = TRUE,                # With downstream clustering of reads
    min_muts = min_muts,        # Skip karyotypes with less then these muts
    QC_type = "D",              # QC with the deconvolution classifier
    N_max = N_max,              # Downsample a karyotype if too many muts
    ...
  )

  # What has not been fit
  which_null = sapply(
    all_fits,
    function(x) all(is.null(x$mobster)) & all(is.null(x$bmix))
  )

  all_fits = all_fits[!which_null]
  subset_mobster_fits = lapply(all_fits, function(x) x$mobster)
  subset_bmix_fits = lapply(all_fits, function(x) x$bmix)

  # 3) Results assembly
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()
  results$type = "Deconvolution pipeline with raw VAF and karyotypes"
  class(results) = "evopipe_rawk"

  # Special case ~ nothing to time, annotate it
  if(length(all_fits) > 0) results$with_fits = TRUE
  else results$with_fits = FALSE

  # Table clustering assignments
  assignment_table = lapply(names(all_fits), function(k)
  {
    x = all_fits[[k]]
    evoverse:::deconvolution_table_assignments(M = x$mobster, B = x$bmix) %>%
      dplyr::mutate(karyotype = k)
  })
  assignment_table = Reduce(dplyr::bind_rows, assignment_table)

  # Table summary fits and MOBSTER QC
  summary_table = lapply(names(all_fits), function(k)
  {
    x = all_fits[[k]]
    evoverse:::deconvolution_table_summary(M = x$mobster, B = x$bmix) %>%
      dplyr::mutate(karyotype = k)
  })
  summary_table = Reduce(dplyr::bind_rows, summary_table)

  # Complete the S3 object with fits and input
  results$mobster = subset_mobster_fits
  results$bmix = subset_bmix_fits

  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input
  )

  # Tables
  results$table$summary = summary_table
  results$table$clustering_assignments = assignment_table

  # Determine karyotypes with sufficiently good quality data
  # CNAqc_input = CNAqc::analyze_peaks(CNAqc_input, matching_epsilon = matching_epsilon_peaks)

  # 2) MOBSTER analysis of karyotypes ~ all analysed
  # mobster_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
  #   mutations = CNAqc_input$snvs,
  #   karyotypes = karyotypes,
  #   min_muts = min_muts,
  #   QC_type = "D",
  #   ...
  # )

  # 3) BMix -- Binomial model -- on each karyotype, non-tail mutations
  # bmix_best_fits = evoverse:::deconvolution_nonread_counts_bmix_VAF(x = mobster_fits,
  #                                                                   karyotypes = karyotypes,
  #                                                                   min_muts =  min_muts)
  # names(bmix_best_fits) = names(mobster_fits)

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

  # # Plot and tables
  # results$table$summary = evoverse:::deconvolution_table_summary(mobster_fits, bmix_best_fits)
  #
  # # Clustering assignments
  # results$table$clustering_assignments = evoverse:::deconvolution_table_assignments(mobster_fits, bmix_best_fits)

  # Data id
  results$description = sample
  results$log = paste0(
    "",
    Sys.time(),
    '. evoverse pipeline for subclonal deconvolution from VAF data and karyotype: with fits ', results$with_fits
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



#' S3 print for the deconvolution pipeline from raw VAF with karyotypes
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.evopipe_rawk = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_ctime'))

  # Print pipeline objects
  cli::cli_rule(paste(crayon::bgYellow(
    crayon::black("[ Evoverse ] {.value {x$description}}")
  ),
  '{.value {x$type}}'))
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$input$cnaqc)

  # PASS/FAIL
#   pass = x$table$summary %>% dplyr::filter(QC == "PASS") %>% dplyr::select(karyotype, starts_with('QC'))
#   fail = x$table$summary %>% dplyr::filter(QC == "FAIL") %>% dplyr::select(karyotype, starts_with('QC'))
#
#   if (nrow(pass) > 0) {
#     cli::cli_rule(crayon::bgGreen(" QC PASS "),
#                   right = paste0("PASS rate (%): ", x$QC$f_PASS))
#     print(pass)
#   }
#
#   if (nrow(fail) > 0) {
#     cli::cli_rule(crayon::bgRed(" QC PASS "),
#                   right = paste0("FAIL rate (%): ", 100 - x$QC$f_PASS))
#     print(fail)
#   }
}

#' S3 plot for the deconvolution pipeline from raw VAF with karyotypes
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.evopipe_rawk = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_qc'))

  # Figure assembly
  # mobster_fits = x$mobster
  # bmix_fits = x$bmix
  # cna_obj = x$input$cnaqc
  #
  # qc_table = x$table$summary
  #
  # # 2 x 1 CNAqc plot
  # groups = names(mobster_fits)
  # empty_panel = CNAqc:::eplot()
  #
  # # MOBSTER plots, sourrounded by a coloured box by QC
  # mob_fits_plot = lapply(groups,
  #                        function(y)
  #                        {
  #                          if (all(is.null(mobster_fits[[y]])))
  #                            return(empty_panel)
  #                          return(CNAqc:::eplot())
  #
  #                          qc_entries = qc_table %>% dplyr::filter(karyotype == !!y)
  #
  #                          qc = ifelse(qc_entries$QC == "FAIL", "indianred3", 'forestgreen')
  #
  #                          mobster::plot.dbpmm(mobster_fits[[y]]) +
  #                            labs(title = bquote("QC " ~ .(qc_entries$QC) ~ p["PASS"] ~
  #                                                  '=' ~ .(qc_entries$QC_prob))) +
  #                            theme(title = element_text(color = qc),
  #                                  panel.border = element_rect(colour = qc,
  #                                                              fill = NA))
  #                        })
  #
  #
  # # CNA plot
  # cna_plot = empty_panel
  # if (!is.null(cna_obj))
  #   cna_plot = CNAqc::plot_segments(cna_obj, circular = TRUE)
  #
  # # Top panel: CNA + MOBSTER
  # mob_fits_plot =  append(list(cna_plot), mob_fits_plot)
  # figure = ggpubr::ggarrange(
  #   plotlist = mob_fits_plot,
  #   nrow = 1,
  #   ncol = length(groups) + 1,
  #   labels = c("CNA", groups)
  # )
  #
  # # Set the figure title and captions
  # figure = ggpubr::annotate_figure(
  #   figure,
  #   top = ggpubr::text_grob(
  #     bquote(bold("Dataset. ") ~ .(x$description)),
  #     hjust = 0,
  #     x = 0,
  #     size = 15
  #   ),
  #   bottom = ggpubr::text_grob(
  #     bquote(.(x$log)),
  #     hjust = 0,
  #     x = 0,
  #     size = 10
  #   )
  # )
  #
  # return(figure)
}
