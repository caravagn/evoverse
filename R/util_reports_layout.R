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

# devtools::install_github("tarkomatas/myprettyreport")

# Raw data panel
get_raw_data_panel = function(x)
{
  USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')

  # Raw mutation data
  raw_muts = x$snvs %>% dplyr::mutate(karyotype = ifelse(karyotype %in% USE_KARYOTYPES, karyotype, "other"))
  ccf_data = Reduce(dplyr::bind_rows,
                    lapply(x$CCF_estimates, function(x)
                      x$mutations)) %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% USE_KARYOTYPES, karyotype, "other"))

  karyo_colors = CNAqc:::get_karyotypes_colors(c(USE_KARYOTYPES, 'other'))

  cvgp = ggplot(data = raw_muts, aes(DP, fill = karyotype)) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Sequencing depth (coverage)",
         caption = paste0("Median coverage ", median(x$snvs$DP), 'x')) +
    geom_vline(
      xintercept = median(raw_muts$DP),
      color = 'orange',
      linetype = 'dashed',
      size = .3
    ) +
    scale_fill_manual(values = karyo_colors)

  vafp = ggplot(data = raw_muts, aes(VAF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, 1) +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Raw VAF",
         caption = paste0("VAF < 0.05 (5%) = ", sum(x$snvs$VAF < 0.05))) +
    scale_fill_manual(values = karyo_colors)

  ccfp = ggplot(data = ccf_data,
                aes(CCF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, max(ccf_data$CCF, na.rm = T) %>% ceiling) +
    CNAqc:::my_ggplot_theme() +
    labs(title = "CCF values",
         caption = paste0("CCF < 0.05 (5%) = ", sum(ccf_data$CCF < 0.05, na.rm = T))) +
    scale_fill_manual(values = karyo_colors)

  return(list(cvgp, vafp, ccfp))
}

# Annotate title on top of a plot
annot_t = function(x, t) {
  ggpubr::annotate_figure(x,
                          top = ggpubr::text_grob(
                            t,
                            hjust = 0,
                            x = 0,
                            face = 'plain'
                          ))

}

# Add margins to a plot
fun_marg = function(x, fact, w = 1, h = 1) {
  x + theme(plot.margin = unit(c(1 * h, 1 * w, 1 * h, 1 * w) * fact, "cm"))
}


page_one = function(x) {
  N = ifelse(x$n_snvs > 3000, 3000, x$n_snvs)

  # Segments
  A = cowplot::plot_grid(
    CNAqc::plot_depth(x, N = N) + labs(
      title = paste0(
        "Ploidy ",
        x$ploidy,
        ', most prevalent karyotype is ',
        CNAqc:::get_prevalent_karyotype(x),
        " (Major:minor).\n"
      )
    ),
    CNAqc::plot_segments(x),
    align = 'v',
    nrow = 2,
    rel_heights = c(.35, .8)
  )

  Bpl  = evoverse:::get_raw_data_panel(x)

  B = Bpl[[1]]
  C = Bpl[[2]]
  D = Bpl[[3]]

  W = ggpubr::ggarrange(C,
                        B,
                        D,
                        ncol = 3,
                        common.legend = T,
                        legend = 'bottom')

  # Genome distribution and QC
  E = ggpubr::ggarrange(
    CNAqc::plot_karyotypes(x),
    CNAqc::plot_karyotypes(x, type = 'number'),
    common.legend = TRUE,
    legend = 'bottom'
  )

  G = CNAqc::plot_qc(x) + labs(title = paste0("Purity/ ploidy (peaks) and CCF estimates (% assignable)"))


  H = cowplot::plot_grid(
    E,
    G,
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("C", 'D'),
    align = 'h',
    axis = 'x'
  )

  ggpubr::ggarrange(
    A,
    ggplot() + geom_blank(),
    W,
    ggplot() + geom_blank(),
    H,
    ncol = 1,
    heights = c(5.8, 1, 4, 1, 4.5),
    labels = c("A", '', "B", '', "C")
  )
}


page_two = function(x)
{
  sc = getOption('CNAqc_cex', default = 1)

  A = cowplot::plot_grid(
    ggplot() + theme_void(base_size = 10  * sc) + labs(
      title = paste0(
        "Purity and ploidy QC via peak detection (purity p = ",
        x$purity * 100,
        "%)\n"
      )
    ),
    CNAqc::plot_peaks_analysis(x),
    align = 'v',
    nrow = 2,
    rel_heights = c(.1, .95)
  )

  N = x$n_snvs * 0.01
  D = median(x$cna$length) * .1

  segment_ids = x$cna %>%
    filter(n > !!N, length > D) %>%
    pull(segment_id)

  if(length(segment_ids) == 0) B = CNAqc:::eplot()
  else
  {
    B = CNAqc::inspect_segment(x, n = N, l = D) +
      labs(
        title = paste0("Segments with >", N, " mutations (1% of total) and >", D, ' bases (10% of median length)\n'),
        subtitle = NULL
      )
  }

  C1 = CNAqc::plot_arm_fragmentation(x) +
    labs(title = paste0("Arm-level segments fragmentation (% relative to arm length)\n"))
  C2 = CNAqc::plot_segment_size_distribution(x) +
    labs(title = paste0("Segments size distribution\n"))

  cowplot::plot_grid(
    A,
    ggplot() + geom_blank(),
    B,
    ggplot() + geom_blank(),
    ggpubr::ggarrange(C1, C2),
    ncol = 1,
    rel_heights = c(4, 1, 4, 1, 4),
    labels = c("A", '', "B", '', "C"),
    scale = c(.9, 1, 1, 1, 1),
    align = 'v',
    axis = 'y'
  )

}

page_three = function(x)
{
  CNAqc::plot_CCF(x) %>%
    annot_t(t = paste0(
      "QC PASS (per karyotype) if we can assign >90% of mutation burden\n"
    ))
}

# Report multi-page style for the Data QC pipeline
report_multipage_cnaqc_pipeline = function(x, f, cex = .7, sample, score, collate = F)
{
  # Set global cex value for all CNAqc plots
  options(CNAqc_cex = cex)
  require(CNAqc)
  require(patchwork)

  # QC
  PASS = score
  l_PASS = cut(PASS, seq(0, 100, by = 20)) %>% as.numeric
  qc_col = brewer.pal('Spectral', n = 5)[l_PASS]

  # Page one -- segments | vaf | qc
  first_page = evoverse:::page_one(x) %>% evoverse:::fun_marg(fact = 1, w = 1, h = )

  # Page two -- Peak analysis and fragmentation
  second_page = evoverse:::page_two(x) %>% evoverse:::fun_marg(fact = 1, w = 1, h = 2)

  # Page two -- Peak analysis and fragmentation
  third_page = evoverse:::page_three(x) %>% evoverse:::fun_marg(fact = 1, w = 1, h = 2)

  # Report assembly
  require(myprettyreport)
  require(magick)

  # tool <-
  #   image_read("https://caravagn.github.io/CNAqc/reference/figures/logo.png")
  tool = NA

  sample_logo = tool
  # frink <- image_read("https://jeroen.github.io/images/frink.png")
  # sample_logo = c(tool, frink)

  cp_e_params = function() {
    B = brewer.pal('Spectral', n = 5)


    # for(i in seq_along(B))
    # {
    #   grid::grid.rect(gp = grid::gpar(fill = B[i]),
    #                   vp = grid::viewport(layout.pos.row = 2,
    #                                       layout.pos.col = 2),
    #                   x = i,
    #                   y = 2,
    #                   width = 0.2, height = 0.2)
    # }
#
# #
    grid::grid.rect(
      gp = grid::gpar(fill = alpha(qc_col, .3)),
      vp = grid::viewport(layout.pos.row = 2,
                          layout.pos.col = 1),
      width = 0.8,
      height = 0.2
    )
#
    grid::grid.text(paste0("PASS rate ", score), gp=grid::gpar(col=qc_col),
                    vp = grid::viewport(layout.pos.row = 2,
                                        layout.pos.col = 1))
  }

  p1_tit = 'Data quality check (QC) from matched Copy Number and mutation data.'
  p2_tit = 'Purity/ ploidy QC (A), data (B) and fragmentation (C).'
  p3_tit = 'Cancer Cell Fractions estimation and QC.'

  # Special report - 1 page wide collated
  if(collate == FALSE)
  {
    message("Saving single-page collated PDF report")

    ggpubr::ggarrange(
      first_page %>% evoverse:::annot_t(p1_tit),
      second_page %>% evoverse:::annot_t(p2_tit),
      third_page %>% evoverse:::annot_t(p3_tit),
      ncol = 3
    ) %>%
      ggsave(filename = f, height = 250, width = 2.1 * 297, units = 'mm')
  }
  else
  {
    message("Saving multi-page PDF report")

    start_report(size = 'a4', filename = f) %>%
      add_cover_page(
        logo = sample_logo,
        title = sample,
        subtitle = p1_tit,
        logo_size = 0.15,
        theme = 'basic',
        creaton_time = paste0(
          'evoverse pipeline (caravagn.github.io/evoverse), ',
          Sys.time()
        ),
        need_footer = T,
        extra_layout_params = cp_e_params
      ) %>%
      add_new_page(
        plot = first_page,
        title = sample,
        subtitle = 'CNA Segments (A), mutations (B), genome coverage (C) and QC (D).',
        need_header = TRUE,
        header_color = "#f44242",
        logo = sample_logo,
        logo_size = .35,
        theme = "flashy",
        need_footer = TRUE,
        footer_text = "Copyright © Giulio Caravagna | 2020 | evoverse (caravagn.github.io/evoverse)",
        footer_color = "#3d3c3c"
      ) %>%
      add_new_page(
        plot = second_page,
        title = sample,
        subtitle = p2_tit,
        need_header = TRUE,
        header_color = "#f44242",
        logo = sample_logo,
        logo_size = .35,
        theme = "flashy",
        need_footer = TRUE,
        footer_text = "Copyright © Giulio Caravagna | 2020 | evoverse (caravagn.github.io/evoverse)",
        footer_color = "#3d3c3c"
      ) %>%
      add_new_page(
        plot = third_page,
        title = sample,
        subtitle = p3_tit,
        need_header = TRUE,
        header_color = "#f44242",
        logo = sample_logo,
        logo_size = .35,
        theme = "flashy",
        need_footer = TRUE,
        footer_text = "Copyright © Giulio Caravagna | 2020 | evoverse (caravagn.github.io/evoverse)",
        footer_color = "#3d3c3c"
      ) %>%
      end_report()
  }

}


# report_multipage_cnaqc_pipeline(x)
#
# library(pdftools)
#
# pdf_combine(c("report_output.pdf", "report_output.pdf"), output = "joined.pdf")
# pdf_compress('joined.pdf', output = "joined3.pdf")
#
#
#
# plotflow:::mergePDF(
#   in.file = paste(c("report_output.pdf", "report_output.pdf"), collapse=" "),
#   file="joined3.pdf"
# )
#
