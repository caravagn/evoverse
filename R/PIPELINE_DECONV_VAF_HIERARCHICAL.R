#' Title
#'
#' @param x
#' @param karyotypes
#' @param description
#' @param min_muts
#' @param min_vaf
#' @param enforce_QC_PASS
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
pipeline_subclonal_deconvolution_hierarchical_VAF = function(
  x,
  karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
  description = "VAF by karyotype hierarchical deconvolution sample",
  min_muts = 150,
  min_vaf = 0.,
  enforce_QC_PASS = TRUE,
  ...) {

  pio::pioHdr("Evoverse", crayon::italic('Raw VAF by karyotype hierarchical subclonal deconvolution pipeline'))
  cat('\n')

  if(!inherits(x, "evopipe_qc")) stop("Input 'x' should be the output of the evoverse data QC pipeline. See ?pipeline_qc_copynumbercalls.")

  # Load input data
  cli::cli_h1("Input data for sample {.field {description}}")
  cat("\n")

  CNAqc_input = x$cnaqc

  print(x)

  results <-  list()


  all_fits <- mobster::mobsterh_fit(x = x,
                                    n_t = min_muts,
                                    vaf_filter = min_vaf,
                                    karyotypes = karyotypes,
                                    ...)


  if(is.null(all_fits)) results$with_fits = FALSE
  else results$with_fits = TRUE

  results$input <-  CNAqc_input
  results$type = "Deconvolution hierarchical pipeline with raw VAF"
  results$description <- description
  class(results) = "evopipe_rawkh"

  if(results$with_fits) {

    # 3) Results assembly
    cat("\n")
    cli::cli_process_start("Pipeline results assembly")
    cat("\n")


    # Complete the S3 object with fits and input
    results$mobster <-  all_fits$best

    results$table$summary <- list(beta = mobster:::get_beta(all_fits$best),pareto = mobster:::get_pareto(all_fits$best),
                                  mixture_weights = mobster:::get_mixture_weights(all_fits$best),assignment_probs = mobster:::get_assignment_probs(all_fits$best)
                                  )
    results$table$clustering_assignments <- all_fits$best$data
    results$table$QC_table <- x$QC$QC_table
    results$type = "Deconvolution hierarchical pipeline with raw VAF"
    class(results) = "evopipe_rawkh"

    results$QC <- ifelse(any((results$table$QC_table %>%
                            filter(type == "Peaks", karyotype %in% names(all_fits$best$model_parameters)) %>%
                            pull(karyotype) %>%  unique()) %in% "FAIL"),
                         "FAIL", "PASS")
    results$architecture <- ifelse(all_fits$best$run_parameters$K == 0, "Monoclonal", "Polyclonal")

    cat('\n')
    cli::cli_h3("Tumour architecture: {.field {results$architecture}}.")

    cli::cli_process_done()

  }

  results$log = paste0(
    Sys.time(),
    '. evoverse pipeline for deconvolution from raw VAF and karyotypes: with fits ', results$with_fits
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
print.evopipe_rawkh = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_rawkh'))

  # Print pipeline objects
  cli::cli_rule(paste(crayon::bgYellow(
    crayon::black("[ Evoverse ] {.value {x$description}}")
  ),
  '{.value {x$type}}'))
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$input)

  # QC of the analysis
  cat("\n")
  cli::cli_rule("QC table for the analysis")
  cat("\n")

  pio::pioDisp(x$table$QC_table)
  cat('\n')

  # Summary message
  if (!is.null(x$QC) && x$QC == "PASS")
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
plot.evopipe_rawkh = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_rawkh'))

  # Figure assembly
  mobster_fits = x$mobster
  cna_obj = x$input

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)

  karyotypes_list = c("1:0", "1:1" , "2:0", "2:1", "2:2")

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = mobster:::plot.dbpmmh(mobster_fits)

  figure = ggpubr::ggarrange(
    CNAqc::plot_segments(cna_obj, circular = FALSE, highlight = groups) + labs(title = x$description),
    CNAqc::plot_peaks_analysis(cna_obj, empty_plot = F),
    mob_fits_plot,
    nrow = 3,
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


