#' Plot a 2-dimensional VAF scatterplot.
#'
#' @description This function plots a 2D scatterplot
#' of the input VAF values in sample`s1` versus `s2`.
#' This function can also subsample the number of plot
#' points, which helps to allow for a faster rendering
#' of the returned plot. The number of mutations shown
#' is always reported, as well as the proportion with
#' respect to the total number of mutations.
#'
#' @param x An `evoverse` object.
#' @param s1 The first sample name, by default `x$samples[1]` the first sample in the data.
#' @param s1 The second sample name, by default `x$samples[1]` the second sample in the data.
#' @param N Maximum number of points to plot, the overall percentage is reported.
#'
#' @return A `ggplot` object
#'
#' @seealso Function \link{\code{plot_2D_VAF_MOBSTER}} extends this plot visualizing also
#' MOBSTER clusters' results.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Plot the 1st vs the 2nd sample (default parameters)
#' plot_2D_VAF(example_evoverse)
plot_2D_VAF = function(x,
                       s1 = x$samples[1],
                       s2 = x$samples[2],
                       N = 10000)
{
  px = VAF(x, samples = s1)
  py = VAF(x, samples = s2)

  points = full_join(px, py, by = 'id') %>%
    filter(value.x > 0 & value.y > 0)

  N_all = nrow(points)
  if(nrow(points) > N)
  {
    message("N = ", N, ' - using only a subset of the data points.')
    points = points %>% sample_n(N)
  }
  else N = nrow(points)

  label = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  ggplot(points,
         aes(x = value.x, y = value.y)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .31) +
    scale_fill_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
    geom_point(alpha = .5, size = 1) +
    guides(fill = FALSE) +
    labs(title = paste0(s1, ' vs ', s2),
         x = s1,
         y = s2) +
    my_ggplot_theme() +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_rug() +
    annotate("label", fill = 'white', x = .9, y = .9, label = label, size = 2, hjust = 1)

}
