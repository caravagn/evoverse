#' Plot a dataset (complex figure).
#'
#' @description Create a complex figure via \code{ggpubr} assembling several plots. The figure
#' reports 1) the adjusted VAF (full spectrum), 2) the zoomed version of the adjusted
#' VAF for values below 0.20 (low-frequency spectrum), 3) the depth against the VAF
#' of each mutation and 4) the depth distribution.
#'
#'
#' @param x A dataset of class \code{"mbst_data"}.
#' @param scales Scales for ggplot facets, default is \code{'fixed'} but \code{'free'} can also help.
#' @param cex Cex of the plots.
#' @param ...
#'
#' @return A ggpubr figure object.
#' @export
#'
#' @examples
#' TODO
plot.mbst_data = function(x, scales = 'fixed', cex = 1, ...)
{
  stopifnot(inherits(x, "mbst_data"))

  # 4 plots
  pl_Full = plot_VAF(x,
                     scales = scales,
                     cex = cex,
                     range = c(0, Inf))

  pl_LowFreq = plot_VAF(x,
                        scales = scales,
                        cex = cex,
                        range = c(0, 0.2))

  pl_DP = plot_DP(x,
                  scales = scales,
                  cex = cex,
                  range = c(0, Inf))

  pl_DP_V = plot_DP_VAF(x, scales = scales)

  # 1 figure
  figure = ggpubr::ggarrange(pl_Full,
                             pl_LowFreq,
                             pl_DP_V,
                             pl_DP,
                             ncol = 1,
                             nrow = 4)

  figure = ggpubr::annotate_figure(figure, top = x$description)

  return(figure)

}

# =-=-=-=-=-=-=-=-
# Aux methods
# =-=-=-=-=-=-=-=-
plot_VAF = function(x,
                    scales = 'fixed',
                    cex = 1,
                    range = c(0, Inf))
{
  stopifnot(inherits(x, "mbst_data"))

  vaf = VAF(x) %>%
    filter(value > range[1], value < range[2])

  max.VAF = max(vaf$value)
  if (max.VAF < 1 & range[2] > 1)
    max.VAF = 1

  ggplot(vaf, aes(value, fill = sample)) +
    facet_wrap( ~ sample, nrow = 1, scales = scales) +
    theme_light(base_size = 8) +
    guides(fill = FALSE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    geom_vline(aes(xintercept = 0.5), colour = 'black', size = .3) +
    geom_vline(
      aes(xintercept = 0.25),
      colour = 'darkblue',
      linetype = "longdash",
      size = .5
    ) +
    geom_vline(
      aes(xintercept = 0.33),
      colour = 'steelblue',
      linetype = "longdash",
      size = .5
    ) +
    geom_vline(
      aes(xintercept = 0.66),
      colour = 'blue',
      linetype = "longdash",
      size = .3
    ) +
    geom_vline(
      aes(xintercept = 1),
      colour = 'black',
      linetype = "longdash",
      size = .2
    ) +
    geom_histogram(binwidth = 0.01, alpha = .8) +
    # geom_density() +
    xlim(0, max.VAF) +
    labs(
      title = bquote(bold("Adjusted VAF")),
      subtitle = paste0('Plot range: ', range[1], '-', range[2]),
      x = 'Adjusted VAF',
      y = 'Counts'
    )
}

plot_DP = function(x,
                   scales = 'fixed',
                   cex = 1,
                   range = c(0, Inf))
{
  stopifnot(inherits(x, "mbst_data"))

  dp = DP(x) %>% filter(value > range[1], value < range[2])

  ggplot(dp, aes(value, fill = sample)) +
    geom_histogram(binwidth = 5, alpha = .8) +
    facet_wrap( ~ sample, nrow = 1, scales = scales) +
    theme_light(base_size = 8) +
    guides(fill = FALSE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    labs(
      title = bquote(bold("Depth (DP)")),
      subtitle = paste0('Plot range: ', range[1], '-', range[2]),
      x = 'DP',
      y = 'Counts'
    )

}

plot_DP_VAF = function(x, scales = 'fixed', cex = 1)
{
  VAF_val = VAF(x) %>% filter(value > 0)
  DP_val = DP(x) %>% filter(value > 0)

  ggplot(bind_rows(DP_val, VAF_val) %>% spread(variable, value),
         aes(x = DP, y = VAF, color = sample)) +
    geom_point(alpha = 1, size = .4) +
    facet_wrap(~ sample, nrow = 1,  scales = 'free') +
    guides(fill = FALSE, color = FALSE) +
    labs(title = bquote(bold("Depth (DP) versus Adjusted VAF")),
         x = 'Depth (DP)',
         y = 'Adjusted VAF') +
    theme_light(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    scale_x_log10() +
    stat_density2d(alpha = .3)
}
