#' Pipeline to time aneuploidy.
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
                                      min_muts = 150,
                                      min_VAF = 0.05,
                                      description = "Chromosomal timing sample",
                                      N_max = 15000,
                                      ...
)
{
  pio::pioHdr("Evoverse", crayon::italic('Chromosomal timing pipeline'))
  cat('\n')

  #
  # 1) Load input data -- this is a common function to all deconvolution-based pipelines
  #
  cat("\n")
  cli::cli_process_start("Loading input data")
  cat("\n")

  CNAqc_input = evoverse:::deconvolution_prepare_input(mutations, cna, purity, min_VAF = min_VAF)

  print(CNAqc_input)
  cli::cli_process_done()

  #
  # 2) MOBSTER analysis of karyotypes, this returns only the best fit
  #
  mobster_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
    mutations = CNAqc_input$snvs,
    karyotypes = timeable,
    min_muts = min_muts,
    QC_type = "T",
    N_max = N_max,
    ...
  )

  #
  # 3) Results assembly
  #
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  results = list()
  results$type = "Chromosome timing pipeline"
  class(results) = "evopipe_ctime"

  # Fits
  results$mobster = mobster_fits

  # Input
  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input)

  # Plot and tables
  results$table$summary = evoverse:::deconvolution_table_summary(mobster_fits, bmix_fits = NULL)

  # Clustering assignments
  results$table$clustering_assignments = evoverse:::deconvolution_table_assignments(mobster_fits, bmix_fits = NULL)

  # Data id
  results$description = description
  results$log = paste0("", Sys.time(), '. evoverse pipeline for chromosome timing. QC: ', results$table$summary$QC_type[1])

  cli::cli_process_done()


  return(results)
}

#' S3 print for the chromosome timing pipeline results.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.evopipe_ctime = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_ctime'))

  # Print pipeline objects
  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ Evoverse ] {.value {x$description}}")),
      '{.value {x$type}}'
    )
  )
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$input$cnaqc)

  # PASS/FAIL
  pass = x$table$summary %>% dplyr::filter(QC == "PASS") %>% dplyr::select(karyotype, starts_with('QC'))
  fail = x$table$summary %>% dplyr::filter(QC == "FAIL") %>% dplyr::select(karyotype, starts_with('QC'))

  if(nrow(pass) > 0) {
    cli::cli_rule(crayon::bgGreen(" QC PASS "), right = paste0("PASS rate (%): ", x$QC$f_PASS))
    print(pass)
  }

  if(nrow(fail) > 0) {
    cli::cli_rule(crayon::bgRed(" QC PASS "), right = paste0("FAIL rate (%): ", 100 - x$QC$f_PASS))
    print(fail)
  }
}

# S3 plot for the chromosome timing pipeline results.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.evopipe_ctime = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_qc'))

  # Figure assembly
  mobster_fits = x$mobster
  bmix_fits = x$bmix
  cna_obj = x$input$cnaqc

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)
  empty_panel = CNAqc:::eplot()

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(
    groups,
    function(y)
    {
      if (all(is.null(mobster_fits[[y]]))) return(empty_panel)
        return(CNAqc:::eplot())

      qc_entries = qc_table %>% dplyr::filter(karyotype == !!y)

      qc = ifelse(qc_entries$QC == "FAIL", "indianred3", 'forestgreen')

      mobster::plot.dbpmm(mobster_fits[[y]]) +
        labs(title = bquote("QC " ~ .(qc_entries$QC) ~ p["PASS"] ~
                              '=' ~ .(qc_entries$QC_prob))) +
        theme(title = element_text(color = qc),
              panel.border = element_rect(
                colour = qc,
                fill = NA
              ))
    })


  # CNA plot
  cna_plot = empty_panel
  if(!is.null(cna_obj)) cna_plot = CNAqc::plot_segments(cna_obj, circular = TRUE)

  # Top panel: CNA + MOBSTER
  mob_fits_plot =  append(list(cna_plot), mob_fits_plot)
  figure = ggpubr::ggarrange(
    plotlist = mob_fits_plot,
    nrow = 1,
    ncol = length(groups) + 1,
    labels = c("CNA", groups)
  )

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(bquote(bold("Dataset. ") ~ .(x$description)), hjust = 0, x = 0, size = 15),
    bottom = ggpubr::text_grob(bquote(.(x$log)), hjust = 0, x = 0, size = 10)
  )

  return(figure)
}
