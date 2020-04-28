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
#' x = pipeline_subclonal_deconvolution_CCF(mutations, cna = cna, purity = purity,  N_max = 1000, auto_setup = 'FAST')
#' print(x)
pipeline_subclonal_deconvolution_CCF = function(mutations,
                                                cna,
                                                purity,
                                                karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                                ccf_method = "ENTROPY",
                                                cutoff_QC_PASS = 0.25,
                                                reference = 'GRCh38',
                                                min_VAF = 0.05,
                                                min_muts = 50,
                                                description = "Example CCF deconvolution",
                                                N_max = 15000,
                                                ...)
{
  pio::pioHdr("Evoverse",
              crayon::italic('~ Subclonal deconvolution pipeline for CCF data'))
  cat('\n')

  #
  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  cli::cli_h1("Loading input data for sample {.field {description}}")
  cat("\n")

  CNAqc_input = evoverse:::deconvolution_prepare_input(mutations,
                                                       cna,
                                                       purity,
                                                       reference = reference,
                                                       min_VAF = min_VAF)

  print(CNAqc_input)

  # Return object will contain input data
  results = list()
  results$type = "Deconvolution pipeline with CCFs"
  class(results) = "evopipe_ccf"

  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input
  )

  results$description = description

  # Check what is actually available in cna_obj
  which_karyo = CNAqc_input$n_karyotype[karyotypes]
  which_karyo = which_karyo[!is.na(which_karyo) & which_karyo > min_muts]
  which_karyo = names(which_karyo)

  # Force exit if there is no suitable karyotype among the ones required
  if(length(which_karyo) == 0)
  {
    results$with_fits = FALSE
    results$mobster = results$bmix = results$table$clustering_assignments = results$table$summary = NULL

    # Data id
    results$log = paste0(
      Sys.time(),
      '. evoverse pipeline for subclonal deconvolution from CCFs: with fits ', results$with_fits
    )

    return(results)
  }

  cli::boxx(
    paste0('The pipeline will analyse: ', paste0(which_karyo, collapse = ', ')),
    background_col = "blue",
    col = 'white'
  ) %>% cat()


  # Compute CCF, and determine the QC = PASS CCF estimates
  CNAqc_input = CNAqc::compute_CCF(CNAqc_input,
                                   karyotypes = which_karyo,
                                   cutoff_QC_PASS = cutoff_QC_PASS,
                                   method = ccf_method)

  # Store this upgraded version of mapping
  results$input$cnaqc = CNAqc_input

  # These have PASS CCF
  karyotypes_with_good_CCF = CNAqc:::compute_QC_table(CNAqc_input)$QC_table %>%
    dplyr::filter(type == 'CCF', QC == 'PASS') %>%
    pull(karyotype)

  # These not, we remove them from inside CNAqc_input
  which_bad_karyotypes = setdiff(names(CNAqc_input$CCF_estimates), karyotypes_with_good_CCF)

  if(length(which_bad_karyotypes) > 0){
    cli::cli_alert_warning("Karyotype(s) {.field {which_bad_karyotypes}} have CCF that do not pass QC and will not be used; try to use another CCF computation method.")

    CNAqc_input$CCF_estimates = CNAqc_input$CCF_estimates[karyotypes_with_good_CCF]
  }

  # Force exit if there is no suitable karyotype among the ones required
  if(length(karyotypes_with_good_CCF) == 0)
  {
    results$with_fits = FALSE
    results$mobster = results$bmix = results$table$clustering_assignments = results$table$summary = NULL

    # Data id
    results$log = paste0(
      Sys.time(),
      '. evoverse pipeline for subclonal deconvolution from CCFs: with fits ', results$with_fits
    )

    return(results)
  }

  #
  # 2) MOBSTER analysis of CCF, this adjust CCF for the VAF etc.
  #
  CCF_fit = evoverse:::deconvolution_mobster_CCF(
    x = CNAqc_input,
    BMix = TRUE,                # With downstream clustering of reads
    min_muts = min_muts,        # Skip karyotypes with less then these muts
    QC_type = "D",              # QC with the deconvolution classifier
    N_max = N_max,              # Downsample a karyotype if too many muts
    min_VAF = min_VAF,          # Minimum VAF (adjusted for CCF)
    ...
  )

  # What has not been fit
  is_null_mobster = all(is.null(CCF_fit$mobster))
  is_null_bmix = all(is.null(CCF_fit$bmix))

  # 3) Results assembly
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  # Special case ~ nothing to time, annotate it
  if(is_null_mobster & is_null_bmix) results$with_fits = FALSE
  else results$with_fits = TRUE

  # Table clustering assignments
  assignment_table = evoverse:::deconvolution_table_assignments(M = CCF_fit$mobster, B = CCF_fit$bmix) %>%
      dplyr::mutate(karyotype = "CCF")

  # Table summary fits and MOBSTER QC
  summary_table = evoverse:::deconvolution_table_summary(M = CCF_fit$mobster, B = CCF_fit$bmix) %>%
      dplyr::mutate(karyotype = "CCF")

  # Complete the S3 object with fits and input
  results$mobster = CCF_fit$mobster
  results$bmix = CCF_fit$bmix

  # Tables
  results$table$summary = summary_table
  results$table$clustering_assignments = assignment_table

  # 4) Final QC and tumour clonal architecture
  cat("\n")
  cli::cli_h1("Assessing the tumour architecture")
  cat("\n")

  # Get QC results for each one of the CCF we analysed.
  QC_CNAqc = CNAqc:::compute_QC_table(CNAqc_input)$QC_table %>%
    dplyr::filter(karyotype %in% which_karyo, type == 'CCF') %>%
    dplyr::rename(CNAqc_QC = QC, CNAqc_type = type) %>%
    dplyr::mutate(n = CNAqc_input$n_karyotype[karyotype])

  QC_CNAqc$p = QC_CNAqc$n/sum(QC_CNAqc$n)

  # We assess the clonal architecture inferred with MOBSTER and BMix
  QC_mobster_CCF = results$table$summary %>%
    dplyr::mutate(
      architecture =
        case_when(
          karyotype %in% "CCF" & BMix_K_B  > 1 ~ "Polyclonal",
          karyotype %in% "CCF" & BMix_K_B  == 1 ~ "Monoclonal"
        )
    ) %>%
    dplyr::select(karyotype,
                  QC,
                  QC_type,
                  QC_prob,
                  architecture,
                  K_beta,
                  tail,
                  BMix_K_B) %>%
    dplyr::rename(
      K_mobster = K_beta,
      K_BMix = BMix_K_B,
      mobster_QC = QC,
      mobster_QC_prob = QC_prob,
      mobster_QC_type = QC_type
    )

  results$table$QC_table = QC_mobster_CCF
  results$table$QC_CCF = QC_CNAqc

  # Marginalise the architecture, the one in the top row
  # is the selected one (Monoclonal or polyclonal)
  Architecture_table = QC_mobster_CCF %>%
    dplyr::filter(mobster_QC == "PASS")

  results$table$Architecture_table = Architecture_table

  # The final QC result is PASS only if there is at least one row
  # in Architecture_table (which requires having a PASS for both
  # CNAqc and MOBSTER fits).
  results$QC = ifelse(nrow(Architecture_table) > 0, "PASS", "FAIL")
  results$architecture = ifelse(nrow(Architecture_table) > 0, Architecture_table$architecture[1],  NA)

  cat('\n')
  cli::cli_h3("Tumour architecture: {.field {results$architecture}}.")

  cli::cli_process_done()

  results$log = paste0(
    Sys.time(),
    '. evoverse pipeline for deconvolution from CCF: with fits ', results$with_fits
  )

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
print.evopipe_ccf = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_ccf'))

  # Print pipeline objects
  cli::cli_rule(paste(crayon::bgYellow(
    crayon::black("[ Evoverse ] {.value {x$description}}")
  ),
  '{.value {x$type}}'))
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$input$cnaqc)

  # QC of CCFs
  cat("\n")
  cli::cli_rule("QC table for CCFs")
  cat("\n")

  pio::pioDisp(x$table$QC_CCF)
  cat('\n')

  # QC of mobster
  cat("\n")
  cli::cli_rule("QC table for CCFss fit")
  cat("\n")

  pio::pioDisp(x$table$QC_table)
  cat('\n')

  # Summary message
  if(x$QC == "PASS")
    cat(
      crayon::green(clisymbols::symbol$tick),
      "Tumour architecture:",
      case_when(
        x$architecture == "Monoclonal" ~ crayon::blue(x$architecture),
        x$architecture == "Polyclonal" ~ crayon::bgMagenta(x$architecture),
        TRUE ~ "Undetermined"
      ),
      "with",
      crayon::bgGreen(' QC PASS ')
    )
  else
    cat(crayon::red(clisymbols::symbol$cross), "Tumour architecture cannot be assessed", crayon::bgRed(' QC FAIL '))
}

#' S3 plot for the deconvolution pipeline from raw VAF with karyotypes
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.evopipe_ccf = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_ccf'))

  # Figure assembly
  mobster_fits = x$mobster
  bmix_fits = x$bmix
  cna_obj = x$input$cnaqc

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = evoverse:::qc_mobster_plot(mobster_fits)

  # BMix plots
  bmix_fits_plots = CNAqc:::eplot()
  if(!all(is.null(x$bmix)))
  {
    subl = paste0('n = ',
                  nrow(x$bmix$input),
                  ", ",
                  paste(
                    names(x$bmix$pi),
                    ' ',
                    round(x$bmix$pi * 100),
                    '%',
                    collapse = ', ',
                    sep = ''
                  ))

    bmix_fits_plots = BMix::plot_clusters(x$bmix, data = x$bmix$input %>% select(NV, DP)) +
      scale_fill_brewer(palette = 'Set2') +
      labs(title = subl, subtitle = NULL)

  }

  # Figure assembly
  figure = ggpubr::ggarrange(
    ggpubr::ggarrange(
      CNAqc::plot_segments(cna_obj, circular = TRUE, highlight = groups) + labs(title = x$description),
      mob_fits_plot,
      bmix_fits_plots,
      nrow = 1),
    CNAqc::plot_CCF(cna_obj, strip = TRUE),
    nrow = 2,
    ncol = 1
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

