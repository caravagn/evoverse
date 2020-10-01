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
#' Run \code{vignette("a3_Deconv_vaf_karyo", package = "evoverse")} to see
#' the corresponding vignette.
#'
#' @param x An object returned from function \code{pipeline_qc_copynumbercalls},
#' which refers to the QC analysis of somatic mutations, CNAs and tumour purity.
#' @param min_muts Skip analysing karyotypes with less than \code{min_muts} mutations.
#' @param ... Parameters passed \code{mobster_fit} in package \code{mobster}.
#' @param description A string description of this pipeline run, it will appear
#' in some of the produced plots.
#' @param karyotypes Vector of ..
#' @param N_max
#' @param enforce_QC_PASS If \code{TRUE}, the pipeline will only analyze samples
#' that have PASS status from the data analysis QC
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
#' y = pipeline_subclonal_deconvolution_VAF_karyotype(x, auto_setup = 'FAST', N_max = 500)
#'
#' print(y)
pipeline_subclonal_deconvolution_VAF_karyotype = function(
  x,
  karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
  description = "VAF by karyotype deconvolution sample",
  min_muts = 150,
  N_max = 15000,
  enforce_QC_PASS = TRUE,
  ...)
{
  pio::pioHdr("Evoverse", crayon::italic('Raw VAF by karyotype subclonal deconvolution pipeline'))
  cat('\n')

  if(!inherits(x, "evopipe_qc")) stop("Input 'x' should be the output of the evoverse data QC pipeline. See ?pipeline_qc_copynumbercalls.")

  # Load input data
  cli::cli_h1("Input data for sample {.field {description}}")
  cat("\n")

  CNAqc_input = x$cnaqc

  print(x)

  # Return object will contain input data
  results = list()
  results$type = "Deconvolution pipeline with raw VAF and karyotypes"
  class(results) = "evopipe_rawk"

  results$input = CNAqc_input
  results$description = description

  # Determine what segments can be used. Check the inputs against the QC status inside x
  # only if enforce_QC_PASS = TRUE. Otherwise use all of them, in principle.
  Peaks_entries = QC_peaks = NULL

  if(enforce_QC_PASS)
  {
    Peaks_entries = x$QC$QC_table %>% dplyr::filter(type == "Peaks")
    QC_peaks = Peaks_entries %>% dplyr::filter(QC == "PASS") %>% dplyr::pull(karyotype)
  }
  else
  {
    Peaks_entries = x$QC$QC_table %>% filter(type == "Peaks")
    QC_peaks = Peaks_entries %>% dplyr::filter(!is.na(QC)) %>%  dplyr::pull(karyotype)

    cli::cli_alert_warning("enforce_QC_PASS = FALSE, karyotypes will be used regardless of QC status.")
    print(Peaks_entries %>% dplyr::filter(!is.na(QC)))
  }

  which_karyo = intersect(QC_peaks, karyotypes)

  # Handle special cases where we cannot run mobster
  if(length(which_karyo) == 0)
  {
    reason = case_when(
      is.null(x$QC$QC_table) ~ "There are no QC tables for input 'x', rerun the evoverse data QC pipeline",
      (nrow(Peaks_entries) == 0) ~ "There are no 'Peaks' in the QC tables for input 'x', rerun the evoverse data QC pipeline",
      all(Peaks_entries$QC != "PASS") ~ "All peaks in the input data are failed, there's nothing to compute timings with.",
      TRUE ~ paste0("Unknown error - the following karyotypes are PASS: ", paste(QC_peaks, collapse = ', '), '.')
    )

    cat("\n")
    cat(
      cli::boxx(
        paste0("There is nothing to perform deconvolution here! ",reason),
        padding = 1,
        col = 'white',
        float = 'center',
        background_col = "brown")
    )
    cat("\n")

    # Setup a minimum object
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
    paste0('The pipeline will analyse karyotypes: ', paste0(which_karyo, collapse = ', ')),
    background_col = "blue",
    col = 'white'
  )

  all_fits = evoverse:::deconvolution_mobster_karyotypes_VAF(
    x = CNAqc_input,
    karyotypes = which_karyo,   # Required karyotypes
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

  # Special case ~ nothing to use, annotate it
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

  # Tables
  results$table$summary = summary_table
  results$table$clustering_assignments = assignment_table

  # 4) Final QC and tumour clonal architecture
  cat("\n")
  cli::cli_h1("Assessing the tumour architecture")
  cat("\n")

  # Get QC results for each one of the peaks we analysed. By definition
  # the peaks are OK in the CNAqc object; now we can superimpose the QC
  # status from the deconvolution. We want to
  # - take a decision based only on  karyotypes that pass peak detection and mobster QC
  # - include the karyotype size in the weighted decision
  QC_CNAqc = CNAqc:::compute_QC_table(CNAqc_input)$QC_table %>%
    dplyr::filter(karyotype %in% which_karyo, type == 'Peaks') %>%
    dplyr::rename(CNAqc_QC = QC, CNAqc_type = type) %>%
    dplyr::mutate(n = CNAqc_input$n_karyotype[karyotype])

  QC_CNAqc$p = QC_CNAqc$n/sum(QC_CNAqc$n)

  # We assess the clonal architecture inferred with MOBSTER and BMix
  QC_mobster_karyotypes = results$table$summary %>%
    dplyr::mutate(
      architecture =
        case_when(
          karyotype %in% c("1:0", "1:1") & BMix_K  > 1 ~ "Polyclonal",
          karyotype %in% c("1:0", "1:1") &
            BMix_K  == 1 ~ "Monoclonal",
          karyotype %in% c("2:0", "2:1", "2:2") &
            BMix_K  > 2 ~ "Polyclonal",
          karyotype %in% c("2:0", "2:1", "2:2") &
            BMix_K  <= 2 ~ "Monoclonal"
        )
    ) %>%
    dplyr::select(karyotype,
                  QC,
                  QC_type,
                  QC_prob,
                  architecture,
                  K_beta,
                  tail,
                  BMix_K) %>%
    dplyr::rename(
      K_mobster = K_beta,
      K_BMix = BMix_K,
      mobster_QC = QC,
      mobster_QC_prob = QC_prob,
      mobster_QC_type = QC_type
      )

  # Merge together both CNAqc and evoverse QC tables
  QC_table = dplyr::full_join(QC_CNAqc, QC_mobster_karyotypes, by = 'karyotype') %>%
    dplyr::select(
      karyotype, n, p,
      K_mobster, tail, K_BMix,
      ends_with('QC'),
      architecture,
      dplyr::everything()
    )

  results$table$QC_table = QC_table

  # Marginalise the architecture, the one in the top row
  # is the selected one (Monoclonal or polyclonal)
  Architecture_table = QC_table %>%
    # dplyr::filter(CNAqc_QC == "PASS", mobster_QC == "PASS") %>% Wrong (disregard QC as it is enforced above)
    dplyr::filter(mobster_QC == "PASS") %>% # Correct, what has been done is done..
    dplyr::group_by(architecture) %>%
    dplyr::summarise(p = sum(p), n = sum(n)) %>%
    dplyr::arrange(desc(p))

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
print.evopipe_rawk = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_rawk'))

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
plot.evopipe_rawk = function(x, ...)
{
  stopifnot(inherits(x, 'evopipe_rawk'))

  # Figure assembly
  mobster_fits = x$mobster
  bmix_fits = x$bmix
  cna_obj = x$input

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)

  karyotypes_list = c("1:0", "1:1" , "2:0", "2:1", "2:2")

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(karyotypes_list,
                         function(x)
                         {
                           if (x %in% groups)
                             return(evoverse:::qc_mobster_plot(mobster_fits[[x]]))
                           else
                             return(CNAqc:::eplot())
                         })

  # BMix plots
  bmix_fits_plots = lapply(karyotypes_list,
                         function(x)
                         {
                           if (x %in% groups)
                           {
                             if (all(is.null(bmix_fits[[x]]))) return(CNAqc:::eplot())
                             else{
                               x = bmix_fits[[x]]

                               subl = paste0('n = ',
                                             nrow(x$input),
                                             ", ",
                                             paste(
                                               names(x$pi),
                                               ' ',
                                               round(x$pi * 100),
                                               '%',
                                               collapse = ', ',
                                               sep = ''
                                             ))

                               BMix::plot_clusters(x, data = x$input %>% select(NV, DP)) +
                                 scale_fill_brewer(palette = 'Set2') +
                                 labs(title = subl, subtitle = NULL)
                             }
                           }
                           else
                             return(CNAqc:::eplot())
                         })

  figure = ggpubr::ggarrange(
    CNAqc::plot_segments(cna_obj, circular = FALSE, highlight = groups) + labs(title = x$description),
    CNAqc::plot_peaks_analysis(cna_obj),
    ggpubr::ggarrange(
      plotlist = mob_fits_plot,
      nrow = 1,
      ncol = length(mob_fits_plot),
      labels = names(mob_fits_plot)
    ),
    ggpubr::ggarrange(
      plotlist = bmix_fits_plots,
      nrow = 1,
      ncol = length(bmix_fits_plots),
      labels = names(bmix_fits_plots)
    ),
    nrow = 4,
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
