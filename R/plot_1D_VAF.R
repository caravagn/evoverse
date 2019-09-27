#' Plot a 1-dimensional VAF histogram.
#'
#' @description This function plots a one dimensional
#' histogam of the input VAF, coloured by mutation karyotype,
#' for a required input sample.
#'
#' @param x An `evoverse` object.
#' @param s The sample name, by default `x$samples[1]` is the first sample in the data.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Plot the 1st (default parameters)
#' plot_1D_VAF(example_evoverse)
#'
# Plot the 2nd
#' plot_1D_VAF(example_evoverse, s = example_evoverse$samples[2])
plot_1D_VAF = function(x,
                       s = x$samples[1])
{
  check_is_mobster_mvdata(x)

  ggplot(VAF(x, samples = s) %>% filter(value > 0),
         aes(value, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    labs(title = paste0(s),
         x = 'VAF',
         y = 'n') +
    my_ggplot_theme() +
    xlim(0, 1) +
    geom_rug()
}
