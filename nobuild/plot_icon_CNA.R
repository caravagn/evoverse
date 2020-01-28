#' Title
#'
#' @param x
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
plot_icon_CNA = function(x,
                         colors =
                           c(
                             `0` = 'steelblue',
                             `1` = 'forestgreen',
                             `2` = 'darkorange',
                             `3` = 'red',
                             `3+` = 'darkred'
                           ))
{
  # stopifnot( .....

  # Get CNAqc data for the chromosomes
  segments = lapply(seq_along(x$CNAqc),
                    function(w)
                      CNAqc:::relative_to_absolute_coordinates(x$CNAqc[[w]]$cna) %>%
                      mutate(
                        minor = ifelse(minor > 3, '3+', paste(minor)),
                        Major = ifelse(Major > 3, '3+', paste(Major)),
                        position  = w,
                        sample = names(x$CNAqc)[w]
                      ))

  # Set levels
  segments = Reduce(bind_rows, segments)
  segments$minor = factor(segments$minor, levels = names(colors))
  segments$Major = factor(segments$Major, levels = names(colors))

  # Base plot
  base_plot = CNAqc:::blank_genome() +
    geom_hline(
      yintercept = 1,
      size = .3,
      colour = 'white',
      linetype = 'dashed'
    ) +
    labs(y = 'Sample')


  # if there are 1-1 segments across all samples, shadow them
  one_one = segments %>% group_by(sample) %>%
    filter(Major == '1', minor == '1') %>%
    group_by(chr, from, to) %>%
    filter(n() == length(x$samples)) %>%
    ungroup()

  if (nrow(one_one) > 0)
    base_plot = base_plot +
    geom_rect(
      data = one_one %>% distinct(chr, from, to, .keep_all = TRUE),
      aes(
        xmin = from,
        xmax = to,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .2,
      fill = 'forestgreen'
    )

  # Add major and minor allele
  base_plot = base_plot +
    geom_segment(
      data = segments,
      aes(
        x = from,
        xend = to,
        y = position,
        yend = position,
        color = Major
      ),
      size = 3
    ) +
    geom_segment(
      data = segments,
      aes(
        x = from,
        xend = to,
        y = position - 0.3,
        yend = position - 0.3,
        color = minor
      ),
      size = 3
    ) +
    geom_hline(
      yintercept = 0.3 + seq_along(x$CNAqc),
      size = .3,
      linetype = 'dashed'
    ) +
    scale_color_manual(values = colors) +
    CNAqc:::my_ggplot_theme() +
    scale_y_discrete(limits = x$samples) +
    guides(color = guide_legend("Allele count"))

    cowplot::plot_grid(
      CNAqc::plot_counts(x$CNAqc[[1]]) + labs(subtitle = x$description),
      base_plot,
      nrow = 2,
      ncol = 1,
      align = 'v',
      rel_heights = c(.2, .9)
    )

}
