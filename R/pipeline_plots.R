pipeline_qc_plots = function(x)
{
  segments_plot = CNAqc::plot_segments(x$input)
  # CCF_plot = CNAqc::plot_data_histogram(x$input, which = "CCF")
  peaks_plot = CNAqc::plot_peaks_analysis(x$input)
  ccf_plot = CNAqc::plot_CCF(x$input)

  ggpubr::ggarrange(
    segments_plot,
    # ggpubr::ggarrange(segments_plot, CCF_plot, nrow = 1, widths = c(1, .25)),
    peaks_plot,
    ccf_plot,
    nrow = 3,
    ncol = 1
  )
}
