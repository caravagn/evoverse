#' Plot a 2-dimensional VAF scatterplot, augmented with `MOBSTER` cluster.
#'
#' @description This function plots a 2D scatterplot
#' of the input VAF values in sample`s1` versus `s2`.
#' Compared to function \link{\code{plot_2D_VAF}}, this
#' function colours the points according to `MOBSTER`
#' clusters (with a double label mechanism `A ~ B`
#' meaning that the point is in cluster `A` in one sample
#' and `B` in the other). Similarly to \link{\code{plot_2D_VAF}},
#' this function can also subsample the number of plot
#' points, which helps to allow for a faster rendering
#' of the returned plot. The number of mutations shown
#' is always reported, as well as the proportion with
#' respect to the total number of mutations.
#'
#' @note This function throws an error if the object `x` does
#' not contain `MOBSTER` clusters.
#'
#' @param x A `evoverse` object.
#' @param s1 The first sample name, by default `x$samples[1]` the first sample in the data.
#' @param s1 The second sample name, by default `x$samples[1]` the second sample in the data.
#' @param N Maximum number of points to plot, the overall percentage is reported.
#' @param ... Extra parameters, not used.
#'
#' @return A `ggplot` object plot.
#'
#' @seealso Function \link{\code{plot_2D_VAF}} is a simpler version of this plot
#' which does not visualize `MOBSTER` clusters.
#'
#' @export
#'
#' @examples
#' TODO
plot_2D_VAF_MOBSTER = function(x,
                       s1 = x$samples[1],
                       s2 = x$samples[2],
                       N = 10000,
                       ...)
{
  check_is_mobster_mvdata(x)

  if(!has_mobster_fits(x)) {
    stop("NO MOB")
  }

  points = Clusters(x)

  s1_col = paste0(s1, '.MOBSTER_cluster')
  s2_col = paste0(s2, '.MOBSTER_cluster')

  if(!(s1_col %in% colnames(points)))
    stop("NO MOB", s1)

  if(!(s2_col  %in% colnames(points)))
    stop("NO MOB", s2)

  s1_col_vaf = paste0(s1, '.VAF')
  s2_col_vaf = paste0(s2, '.VAF')

  points = points %>%
    filter(!!s1_col_vaf > 0 & !!s2_col_vaf > 0) %>%
    mutate(
      label = paste0(!!s1_col, ' ~ ', !!s2_col)
    ) %>%
    select(!!s1_col, !!s1_col_vaf, !!s2_col, !!s2_col_vaf, label)

  # Subset if required
  N_all = nrow(points)
  if(nrow(points) > N)
  {
    message("N =", N, ' - using only a subset of the data points.')
    points = points %>% sample_n(N)
  }
  else N = nrow(points)

  label = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  points$label = paste0(
    points %>% pull(!!s1_col), ' ~ ',
    points %>% pull(!!s2_col))

  ggplot(points,
         aes_string(x = s1_col_vaf, y = s2_col_vaf, color = "label")) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .1, color = NA) +
    scale_fill_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
    geom_point(alpha = .5, size = 1) +
    guides(fill = FALSE, color = guide_legend('Cluster')) +
    labs(title = paste0('MOBSTER: ', s1, ' vs ', s2),
         x = s1,
         y = s2) +
    my_ggplot_theme() +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_rug() +
    annotate("label", fill = 'white', x = .9, y = .9, label = label, size = 2, hjust = 1)
}
