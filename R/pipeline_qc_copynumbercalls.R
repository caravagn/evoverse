mutations = CNAqc::example_dataset_CNAqc$snvs
cna = CNAqc::example_dataset_CNAqc$cna
purity = CNAqc::example_dataset_CNAqc$purity
reference = CNAqc::example_dataset_CNAqc$reference

w = pipeline_qc_copynumbercalls(mutations, cna, purity)
w$plot %>% ggsave(filename = 'a.pdf', width = 14, height = 28)

x = w$plot

hd = ggplot_report_header("FAFDF", paragraph = Sys.time())
ww = ggarrange(hd, x, heights = c(0.03, 0.97), ncol = 1, nrow = 2)
ww %>% ggsave(filename = 'a.pdf', width = 14, height = 28)

#' Title
#'
#' @param mutations
#' @param cna
#' @param purity
#' @param reference
#' @param sample
#' @param matching_epsilon_peaks
#'
#' @return
#' @export
#'
#' @examples
pipeline_qc_copynumbercalls = function(
  mutations,
  cna,
  purity,
  reference = 'GRCh38',
  sample = "My Sample",
  matching_epsilon_peaks = 0.025)
{
  pio::pioHdr("Evoverse", italic(paste0('~ Pipeline to determine thq quality of somatic mutation and CNA data')))
  cat('\n')

  x = CNAqc::init(mutations, cna, purity, ref = reference)

  USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
  prevalent_karyo = CNAqc:::get_prevalent_karyotype(x, rank = 1)

  # Segments
  segments = CNAqc::plot_segments(x)
  # segments = ggpubr::ggarrange(
  #   plot_label,
  #   segments,
  #   ncol = 1,
  #   nrow = 2,
  #   heights = c(1, 4)
  # )

  ##### ##### ##### ##### ##### ##### ##### #####
  # peak detection ~ QC purity/ ploidy
  ##### ##### ##### ##### ##### ##### ##### #####
  x = CNAqc::analyze_peaks(x,
                           matching_epsilon = matching_epsilon_peaks,
                           karyotypes = USE_KARYOTYPES)

  adjustment = round(x$peaks_analysis$score * 2, 2)

  peaks_score_title = bquote(
    bold("  Purity/ploidy QC : ") ~ rho ~'= '* .(x$peaks_analysis$score) ~ ' (' * epsilon ~ '=' ~ .(matching_epsilon_peaks) * ')'
  )

  peaks = suppressWarnings(suppressMessages(CNAqc::plot_peaks_analysis(x)))

  ##### ##### ##### ##### ##### ##### ##### #####
  # CCF values
  ##### ##### ##### ##### ##### ##### ##### #####
  x = CNAqc::compute_CCF(x, karyotypes = USE_KARYOTYPES)

  ccf_plot = lapply(
    USE_KARYOTYPES,
    function(k)
    {
      if(!(k %in% names(x$CCF_estimates))) return(CNAqc:::eplot())
      suppressWarnings(CNAqc:::plot_mutation_multiplicity_entropy(x, k))
    }
  )

  ccf_plot = ggpubr::ggarrange(plotlist = ccf_plot, nrow = length(USE_KARYOTYPES), ncol = 1)

  ccf_plot = ggpubr::annotate_figure(ccf_plot,
                                  top = ggpubr::text_grob(x = 0, hjust = 0,  label = "  Cancer Cell Fractions (estimated)")
  )

  ##### ##### ##### ##### ##### ##### ##### #####
  # Fragmentation (after smoothing)
  ##### ##### ##### ##### ##### ##### ##### #####
  x = CNAqc::smooth_segments(x)
  x = CNAqc::detect_arm_overfragmentation(x)
  frags_plot = suppressWarnings(CNAqc::plot_arm_fragmentation(x, zoom = 0))

  x =  CNAqc::detect_wg_overfragmentation(x)


  ##### ##### ##### ##### ##### ##### ##### #####
  # Figure assembly
  ##### ##### ##### ##### ##### ##### ##### #####
  rpass = function(f) {
    ggpubr::annotate_figure(f, top = '')
  }

  figure =
    suppressWarnings(
      ggpubr::ggarrange(
        segments %>% rpass,
        peaks  %>% rpass,
        ccf_plot  %>% rpass,
        frags_plot  %>% rpass,
        heights = c(1, 1, length(USE_KARYOTYPES), 1),
        labels = LETTERS[1:4],
        ncol = 1,
        nrow = 4
      )
    )

  # Report label
  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(x = 0, hjust = 0,
                            label = bquote("  Input purity" ~.(purity) ~ "; CNAqc ploidy " ~ .(x$ploidy) * ' (' * .(prevalent_karyo) * ' karyotype)'),
                            size = 12)
  )

  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(x = 0, hjust = 0,
                            label = bquote(bold("  evoverse CNA/ mutation report : ") ~.(sample)), size = 15)
  )

  # Return
  fit = list()
  fit$fit = x
  fit$plot = figure

  return(fit)
}

