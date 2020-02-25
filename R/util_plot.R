my_ggplot_theme = function(cex = 1)
{
  theme_light(base_size = 10 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}

my_ggplot_karyotype_colors = function()
{
 c(
  `0:0` = 'darkblue',
  `1:0` = 'steelblue',
  `1:1` = ggplot2::alpha('forestgreen', .8),
  `2:0` = 'turquoise4',
  `2:1` = ggplot2::alpha('orange', .8),
  `2:2` = 'indianred3',
  `other` = 'darkred'
  )
}
