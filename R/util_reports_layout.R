ggplot_report_header = function(sample, paragraph = Sys.time())
{
  blank = ggplot() + theme_void()

  button = function(color, radio, size, font = 'plain')
  {
    ggplot(data = data.frame(
      x = 0,
      y = 0,
      label = radio,
      stringsAsFactors = FALSE
    ),
    aes(x = x , y = y, label = label)) +
      theme_void() +
      theme(plot.background = element_rect(fill = color)) +
      geom_text(fontface = font, hjust = 0)
  }

  tline = ggpubr::ggarrange(
    button('gray', sample, size = 10, font = 'bold'),
    blank,
    ncol = 2,
    nrow = 1,
    widths = c(.7, .3)
  )

  dtpar = button('gainsboro', paragraph, size = 10)

  figure = ggpubr::ggarrange(
    blank,
    tline,
    dtpar,
    blank,
    ncol = 1,
    nrow = 4,
    heights = c(.15, 1, .7, .2)
  )
}


devtools::install_github("tarkomatas/myprettyreport")

library(ggplot2)
sample_plot <- ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point() +
  stat_smooth(method = 'lm', color = "#f44242", fill = "#fd9068")

library(magick)
sample_logo <- image_read("https://caravagn.github.io/CNAqc/reference/figures/logo.png")


header_front = paste0(
  'Data quality check (QC) from matched Copy Number and mutation data.'
)



segments_plot = suppressMessages(suppressWarnings(
  cowplot::plot_grid(
    plot_counts(x),
    # plot_vaf(x, N = 10000),
    plot_depth(x, N = 5000),
    plot_segments(x),
    align = 'v',
    nrow = 4,
    rel_heights = c(.15, .15, .8)
  )
))

first_page = suppressMessages(suppressWarnings(
  cowplot::plot_grid(
    # button(color = 'white', radio = ""),
    segments_plot + labs(title = "CNA segments") + theme(plot.margin=unit(c(3,1,1,1) * .3,"cm")),
    button(color = 'white', radio = paste0("Peak detection for purity ", x$purity * 100, '%')),
    peaks + theme(plot.margin=unit(c(2,1,2,1) * .3,"cm")),
    # button(color = 'white', radio = ""),
    cowplot::plot_grid(
      plot_arm_fragmentation(x),
      plot_segment_size_distribution(x),
      ncol = 2,
      rel_widths = c(2,1),
      align = 'v',
      axis = 'y'
      ) + theme(plot.margin=unit(c(2,1,3,1) * .3,"cm")),
    nrow = 4,
    ncol = 1,
    rel_heights = c(1.5, .05, 1, 1.5)
    # scale = c(1, .9, 1, .8, 1, .9)
  )))

second_page = suppressMessages(suppressWarnings(
  cowplot::plot_grid(
    ccf_plot,
    nrow = 1,
    ncol = 1
    # scale = c(1, .9, 1, .8, 1, .9)
  )))

#
# library(myprettyreport)
start_report() %>%
  add_cover_page(
    logo = sample_logo,
    title = sample,
    subtitle = header_front,
    logo_size = 0.15,
    need_footer = T
  ) %>%
  add_new_page(
    plot = first_page,
    title = sample,
    subtitle = 'Segments data and purity/ ploidy QC.',
    need_header = TRUE,
    logo = sample_logo,
    logo_size = 0.1
  ) %>%
  add_new_page(
    plot = second_page,
    title = sample,
    subtitle = 'Cancer Cell Fractions data',
    need_header = TRUE,
    logo = sample_logo,
    logo_size = 0.1
  ) %>%
  end_report()


