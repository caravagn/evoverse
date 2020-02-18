# mutations = CNAqc::example_dataset_CNAqc$snvs
# cna = CNAqc::example_dataset_CNAqc$cna
# purity = CNAqc::example_dataset_CNAqc$purity
# reference = CNAqc::example_dataset_CNAqc$reference
#
# w = pipeline_qc_copynumbercalls(mutations, cna, purity)
# w$plot %>% ggsave(filename = 'a.pdf', width = 15, height = 23)

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
  x = CNAqc::analyze_peaks(x, matching_epsilon = matching_epsilon_peaks, min_density = 0.05, karyotypes = USE_KARYOTYPES)

  adjustment = round(x$peaks_analysis$score * 2, 2)

  peaks_score_title = bquote(
    bold("  Purity/ploidy QC : ") ~ rho ~'= '* .(x$peaks_analysis$score) ~ ' (' * epsilon ~ '=' ~ .(matching_epsilon_peaks) * ')'
  )

  color_peaks_qc = 'forestgreen'
  if(adjustment > 0.03) color_peaks_qc = 'darkgoldenrod2'
  if(adjustment > 0.05) color_peaks_qc = 'chocolate3'
  if(adjustment > 0.10) color_peaks_qc = 'brown3'

  peaks = CNAqc::plot_peaks_analysis(x)
  peaks = ggpubr::annotate_figure(peaks,
                            top = ggpubr::text_grob(x = 0, hjust = 0,  label = peaks_score_title, color = color_peaks_qc)
  )

  ##### ##### ##### ##### ##### ##### ##### #####
  # CCF values
  ##### ##### ##### ##### ##### ##### ##### #####
  x = CNAqc::compute_CCF(x, karyotypes = USE_KARYOTYPES)

  ccf_plot = lapply(
    USE_KARYOTYPES,
    function(k)
    {
      if(!(k %in% names(x$CCF_estimates))) return(CNAqc:::eplot())
      CNAqc:::plot_mutation_multiplicity_entropy(x, k)
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
  frags_plot = CNAqc::plot_arm_fragmentation(x, zoom = 0)

  x =  CNAqc::detect_wg_overfragmentation(x)

  ##### ##### ##### ##### ##### ##### ##### #####
  # Figure assembly
  ##### ##### ##### ##### ##### ##### ##### #####
  rpass = function(f) {
    ggpubr::annotate_figure(f, top = '')
  }

  figure =
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

