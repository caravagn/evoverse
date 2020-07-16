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
#' @param x An object result from the evoverse data QC pipeline, which must contain results
#' from peak detection analysis to decide which CNA events can be timed. This object contains
#' also somatic mutations and CNA data; see function \code{pipeline_qc_copynumbercalls}.
#' @param karyotypes Karyotypes that should be timed, expressed as "Major:minor" notation. Default
#' are the following states \code{timeable = c('2:0', '2:1', '2:2')}. This list of karyotypes will
#' be checked against the list of peaks that PASS QC in \code{x}; only karyotypes that pass and are
#' listed will be used.
#' @param min_muts Karyotypes with less than \code{min_muts} mutations will not be analysed, even
#' if they PASS all the data metrics in \code{x} and are in \code{karyotypes}.
#' @param ... Parameters passed to \code{mobster_fit} in package \code{mobster} to perform the
#' signal deconvolution. Each run is atteempted twice if it FAILS at first; the first run uses
#' parameters specified by this ellipsis, the second run uses the default parameters of \code{mobster_fit}.
#'
#' @seealso pipeline_qc_copynumbercalls
#'
#' @return
#' @export
#'
#' @examples
#' # We use data released with the CNAqc package
#' x = pipeline_qc_copynumbercalls(
#'   mutations = CNAqc::example_dataset_CNAqc$snvs,
#'   cna = CNAqc::example_dataset_CNAqc$cna,
#'   purity = CNAqc::example_dataset_CNAqc$purity
#'   )
#'
#' print(x)
#'
#' # We use x to run the pipeline
#' y = pipeline_chromosome_timing(x, auto_setup = 'FAST', N_max = 500)
#'
#' print(y)
pipeline_chromosome_timing =
  function(x,
           karyotypes = c('2:0', '2:1', '2:2'),
           min_muts = 150,
           description = "Chromosomal timing sample",
           N_max = 15000,
           ...)
  {
    pio::pioHdr("Evoverse", crayon::italic('Chromosomal timing pipeline'))
    cat('\n')

    if (!inherits(x, "evopipe_qc"))
      stop(
        "Input 'x' should be the output of the evoverse data QC pipeline. See ?pipeline_qc_copynumbercalls."
      )

    # Load input data
    cli::cli_h1("Input data for sample {.field {description}}")
    cat("\n")

    CNAqc_input = x$cnaqc

    print(x)

    # Return object will contain input data
    results = list()
    results$type = "Chromosome timing pipeline"
    class(results) = "evopipe_ctime"

    results$input = CNAqc_input
    results$description = description

    # Determine what segments can be used. Check the required inputs against the QC status inside x.
    Peaks_entries = x$QC$QC_table %>% filter(type == "Peaks")
    QC_peaks = Peaks_entries %>% filter(QC == "PASS") %>% pull(karyotype)

    which_karyo = intersect(QC_peaks, karyotypes)

    # Handle special cases where we cannot time
    if (length(which_karyo) == 0)
    {
      reason = case_when(
        is.null(x$QC$QC_table) ~ "There are no QC tables for input 'x', rerun the evoverse data QC pipeline",
        (nrow(Peaks_entries) == 0) ~ "There are no 'Peaks' in the QC tables for input 'x', rerun the evoverse data QC pipeline",
        all(Peaks_entries$QC != "PASS") ~ "All peaks in the input data are failed, there's nothing to compute timings with.",
        TRUE ~ paste0(
          "Unknown error - the following karyotypes are PASS: ",
          paste(QC_peaks, collapse = ', '),
          '.'
        )
      )

      cat("\n")
      cat(
        cli::boxx(
          paste0("There is nothing to time here! ", reason),
          padding = 1,
          col = 'white',
          float = 'center',
          background_col = "brown"
        )
      )
      cat("\n")

      # Setup a minimum object
      results$with_fits = FALSE
      results$mobster = results$bmix = results$table$clustering_assignments = results$table$summary = NULL

      # Data id
      results$log = paste0(
        Sys.time(),
        '. evoverse pipeline for subclonal deconvolution from VAF data and karyotype: with fits ',
        results$with_fits
      )

      return(results)
    }

    cli::boxx(
      paste0(
        'The pipeline will analyse karyotypes: ',
        paste0(which_karyo, collapse = ', ')
      ),
      background_col = "blue",
      col = 'white'
    ) %>% cat()

    all_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
      x = CNAqc_input,
      karyotypes = which_karyo,
      # Required karyotypes
      BMix = FALSE,
      # Withouth downstream clustering of reads
      min_muts = min_muts,
      # Skip karyotypes with less then these muts
      QC_type = "T",
      # QC with the timing classifier
      N_max = N_max,
      # Downsample a karyotype if too many muts
      ...
    )

    # What has not been fit
    which_null = sapply(all_fits,
                        function(x)
                          all(is.null(x$mobster)) &
                          all(is.null(x$bmix)))

    all_fits = all_fits[!which_null]
    subset_mobster_fits = lapply(all_fits, function(x)
      x$mobster)

    # 3) Results assembly
    cat("\n")
    cli::cli_process_start("Pipeline results assembly")
    cat("\n")

    # Special case ~ nothing to time, annotate it
    if (length(all_fits) > 0)
      results$with_fits = TRUE
    else
      results$with_fits = FALSE

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

    # Tables
    results$table$summary = summary_table
    results$table$clustering_assignments = assignment_table

    # Data id
    results$description = description
    results$log = paste0(Sys.time(),
                         '. evoverse pipeline for chromosome timing: with fits ',
                         results$with_fits)

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
  CNAqc:::print.cnaqc(x$input)


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
  cna_obj = x$input

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)
  karyotypes_list = c("2:0", "2:1", "2:2")

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(karyotypes_list,
                         function(x)
                         {
                           if (x %in% groups)
                             return(evoverse:::qc_mobster_plot(mobster_fits[[x]]))
                           else
                             return(CNAqc:::eplot())
                         })

  # Figure
  figure = ggpubr::ggarrange(
    CNAqc::plot_segments(cna_obj, circular = FALSE, highlight = groups) + labs(title = x$description),
    CNAqc::plot_peaks_analysis(cna_obj),
    ggpubr::ggarrange(
      plotlist = mob_fits_plot,
      nrow = 1,
      ncol = length(mob_fits_plot),
      labels = names(mob_fits_plot)
    ),
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

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(figure,
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
                                   ))

  return(figure)
}
