#' Plot a report of the data available for a sample.
#'
#' @description This function creates a multi-plot figure from
#' different visualizations available in `evoverse`, and in the
#' other packages that are `evoverse` interfaces to. The plot
#' regards only one sample in the dataset and shows:
#'
#' \itemize{
#' \item The histogram of the alelleic frequency (VAF) for the sample;
#' \item The histogram of the depth (DP) for the sample;
#' \item The histogram of the number of reads with the variant (NV) for the sample;
#' \item The joint scatter of DP against VAF;
#' \item The copy number segments, plot via the \code{CNAqc} package
#' available at \url{https://github.com/caravagn/CNAqc}.
#' }
#'
#' @param x A mvMOBSTER object.
#' @param sample The sample to plot.
#' @param N_CNA The number of mutations to plot in the CNA segment
#' panel. This parameter is passed to the plotting functions of
#' the \code{CNAqc} package.
#'
#' @return A figure assembled from several \code{ggplot} objects.
#'
#' @import ggpubr
#' @import cowplot
#'
#' @export
#'
#' @examples
#' TODO
plot_sample = function(x, sample, N_CNA = 10000)
{
  # Copy number panel
  CNA = x$CNAqc[[sample]]

  CNA_panel = cowplot::plot_grid(
    plot_counts(CNA),
    plot_vaf(CNA, N = N_CNA),
    plot_depth(CNA, N = N_CNA),
    plot_segments(CNA),
    align = 'v',
    nrow = 4,
    rel_heights = c(.15, .15, .15, .8))

  CNA_panel = annotate_figure(
    CNA_panel,
    top = text_grob("CNA segments", x = 0, h = 0)
  )

  # VAF, DP and NV data
  vp = VAF(x, samples = sample) %>%
    ggplot(aes(value, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, 1) +
    my_ggplot_theme() +
    labs(title = 'VAF')

  dp = DP(x, samples = sample) %>%
    ggplot(aes(value, fill = karyotype)) +
    geom_histogram(bins = 100) +
    my_ggplot_theme() +
    labs(title = 'DP')

  np = NV(x, samples = sample) %>%
    ggplot(aes(value, fill = karyotype)) +
    geom_histogram(bins = 100) +
    my_ggplot_theme() +
    labs(title = 'NV')

  vdp = VAF(x, samples = sample) %>%
    full_join(DP(x, samples = sample), by = 'id') %>%
    ggplot(aes(x = value.x, y = value.y)) +
    geom_point(size = .5, aes(color = karyotype.x)) +
    my_ggplot_theme() +
    guides(color = guide_legend("Karyotype", ncol = 3)) +
    xlim(0, 1) +
    labs(x = 'VAF', y = 'DP') +
    labs(title = 'VAF versus DP')

  mutations_panel = ggarrange(
    vp, dp, np, vdp,
    ncol = 4,
    common.legend = TRUE,
    legend = 'bottom')

  mutations_panel = annotate_figure(
    mutations_panel,
    top = text_grob(paste(sample, ' ~ ', x$description), x = 0, h = 0)
  )

  # Figure
  figure = cowplot::plot_grid(
    mutations_panel,
    CNA_panel,
    align = 'v',
    nrow = 2,
    rel_heights = c(.8, 1)
  )

  figure
}
