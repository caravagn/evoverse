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
#' @param karyotypes Karyotypes that can be timed, expressed as "Major:minor" notation. Default
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
pipeline_chromosome_timing =
  function(mutations,
           cna,
           purity,
           karyotypes = c('2:0', '2:1', '2:2'),
           reference = 'GRCh38',
           min_muts = 150,
           min_VAF = 0.05,
           description = "Chromosomal timing sample",
           N_max = 15000,
           ...)
  {

  pio::pioHdr("Evoverse", crayon::italic('Chromosomal timing pipeline'))
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
  results$type = "Chromosome timing pipeline"
  class(results) = "evopipe_ctime"

  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input
  )

  results$description = description

  #
  # 2) MOBSTER analysis of karyotypes, this returns only the best fit
  #
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
      '. evoverse pipeline for subclonal deconvolution from VAF data and karyotype: with fits ', results$with_fits
    )

    return(results)
  }

  cli::boxx(
    paste0('The pipeline will analyse: ', paste0(which_karyo, collapse = ', ')),
    background_col = "blue",
    col = 'white'
  )

  all_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
    x = CNAqc_input,
    karyotypes = which_karyo,   # Required karyotypes
    BMix = FALSE,               # Withouth downstream clustering of reads
    min_muts = min_muts,        # Skip karyotypes with less then these muts
    QC_type = "T",              # QC with the timing classifier
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

  # 3) Results assembly
  cat("\n")
  cli::cli_process_start("Pipeline results assembly")
  cat("\n")

  # Special case ~ nothing to time, annotate it
  if(length(all_fits) > 0) results$with_fits = TRUE
  else results$with_fits = FALSE

  # Table clustering assignments (Bmix is NULL by definition)
  assignment_table = lapply(names(all_fits), function(k)
  {
    x = all_fits[[k]]
    evoverse:::deconvolution_table_assignments(M = x$mobster, B = NULL) %>%
      dplyr::mutate(karyotype = k)
  })
  assignment_table = Reduce(dplyr::bind_rows, assignment_table)

  # Table summary fits and MOBSTER QC (Bmix is NULL by definition)
  summary_table = lapply(names(all_fits), function(k)
  {
    x = all_fits[[k]]
    evoverse:::deconvolution_table_summary(M = x$mobster, B = NULL) %>%
      dplyr::mutate(karyotype = k)
  })
  summary_table = Reduce(dplyr::bind_rows, summary_table)

  # Complete the S3 object with fits and input
  results$mobster = subset_mobster_fits

  results$input = list(
    mutations = mutations,
    cna = cna,
    purity = purity,
    cnaqc = CNAqc_input
  )

  # Tables
  results$table$summary = summary_table
  results$table$clustering_assignments = assignment_table

  # Data id
  results$description = description
  results$log = paste0(
    Sys.time(),
    '. evoverse pipeline for chromosome timing: with fits ', results$with_fits
    )

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
  cli::cli_rule(paste(crayon::bgYellow(
    crayon::black("[ Evoverse ] {.value {x$description}}")
  ),
  '{.value {x$type}}'))
  cat("\n")

  # Skinnier print
  CNAqc:::print.cnaqc(x$input$cnaqc)


  # PASS/FAIL
  pass = x$table$summary %>% dplyr::filter(QC == "PASS") %>% dplyr::select(karyotype, starts_with('QC'))
  fail = x$table$summary %>% dplyr::filter(QC == "FAIL") %>% dplyr::select(karyotype, starts_with('QC'))

  if (nrow(pass) > 0) {
    cat("\n")
    cli::cli_rule(crayon::bgGreen(" QC PASS "),
                  right = paste0("PASS rate (%): ", x$QC$f_PASS))
    print(pass)
  }

  if (nrow(fail) > 0) {
    cat("\n")
    cli::cli_rule(crayon::bgRed(" QC PASS "),
                  right = paste0("FAIL rate (%): ", 100 - x$QC$f_PASS))
    print(fail)
  }

  cat('\n', crayon::bgBlue(" LOG "), x$log)
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
  stopifnot(inherits(x, 'evopipe_ctime'))

  # Figure assembly
  mobster_fits = x$mobster
  cna_obj = x$input$cnaqc

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(mobster_fits, evoverse:::qc_mobster_plot)

  # CNA plot
  cna_plot = CNAqc::plot_segments(cna_obj, circular = TRUE, highlight = groups) +
    labs(title = x$description)

  # Top panel: CNA + MOBSTER
  figure = ggpubr::ggarrange(
    plotlist = append(list(cna_plot), mob_fits_plot),
    nrow = 1,
    ncol = length(groups) + 1,
    labels = c("", groups)
  )

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(
    figure,
    # top = ggpubr::text_grob(
    #   bquote(bold("Dataset. ") ~ .(x$description)),
    #   hjust = 0,
    #   x = 0,
    #   size = 15
    # ),
    bottom = ggpubr::text_grob(
      bquote(.(x$log)),
      hjust = 0,
      x = 0,
      size = 10
    )
  )

  return(figure)
}
