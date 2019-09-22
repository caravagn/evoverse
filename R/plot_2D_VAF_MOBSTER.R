
#' Plot a 2-dimensional VAF scatterplot, augmented with MOBSTER cluster
#'
#' @description For 2 input samples, this plots a 2D scatterplot
#' of the input VAF in `s1` versus `s2`. The colour of each point
#' represents the result of running MOBSTER on the input data.
#'
#' @param x A `mvMOSTER` object.
#' @param s1 The first sample name, by default `x$samples[1]` the first sample in the data.
#' @param s1 The second sample name, by default `x$samples[1]` the second sample in the data.
#'
#' @return
#' @export
#'
#' @examples
plot_2D_VAF_MOBSTER = function(x,
                       s1 = x$samples[1],
                       s2 = x$samples[2])
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
    mutate(
      label = paste0(!!s1_col, ' ~ ', !!s2_col)
    ) %>%
    select(!!s1_col, !!s1_col_vaf, !!s2_col, !!s2_col_vaf, label)

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
    geom_rug()
}
