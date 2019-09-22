#' Plot a 2-dimensional VAF scatterplot
#'
#' @description For 2 input samples, this plots a 2D scatterplot
#' of the input VAF in `s1` versus `s2`.
#'
#' @param x A `mvMOSTER` object.
#' @param s1 The first sample name, by default `x$samples[1]` the first sample in the data.
#' @param s1 The second sample name, by default `x$samples[1]` the second sample in the data.
#'
#' @return
#' @export
#'
#' @examples
plot_2D_VAF = function(x,
                       s1 = x$samples[1],
                       s2 = x$samples[2])
{
  px = VAF(x, samples = s1)
  py = VAF(x, samples = s2)

  points = full_join(px, py, by = 'id')

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
    geom_rug()
}
