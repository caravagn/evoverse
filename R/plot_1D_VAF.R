#' Plot a 1-dimensional VAF histogram.
#'
#' @description For an input sample, this plots a one dimensional
#' histogam of the input VAF, coloured by mutation karyotype.
#'
#' @param x A `mvMOSTER` object.
#' @param s The sample name, by default `x$samples[1]` the first sample in the data.
#'
#' @return
#' @export
#'
#' @examples
plot_1D_VAF = function(x,
                       s = x$samples[1])
{
  ggplot(VAF(x, samples = s),
         aes(value, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    labs(title = paste0(s),
         x = 'VAF',
         y = 'n') +
    my_ggplot_theme() +
    xlim(0, 1) +
    geom_rug()
}
