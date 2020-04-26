# # Get a QC table
# get_pqc_cna_QC_table = function(x)
# {
#   all_karyptypes = x$peaks_analysis$fits %>% names
#
#   if(all(is.null(x$peaks_analysis)))
#   {
#     peaks_QC = data.frame(
#       karyotype = all_karyptypes,
#       QC = NA,
#       type = 'Peaks',
#       stringsAsFactors = FALSE
#     )
#   }
#   else
#   {
#     peaks_QC = x$peaks_analysis$matches %>%
#       dplyr::select(karyotype, QC) %>%
#       dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
#       dplyr::distinct(karyotype, QC, .keep_all = T) %>%
#       dplyr::arrange(karyotype) %>%
#       dplyr::mutate(
#         value = 1,
#         lab.ypos = cumsum(value) - 0.5 * value,
#         QC = paste(QC),
#         label = karyotype,
#         type = 'Peaks')
#   }
#
#   if(all(is.null(x$CCF_estimates)))
#   {
#     CCF_QC = data.frame(
#       karyotype = all_karyptypes,
#       QC = NA,
#       type = 'CCF',
#       stringsAsFactors = FALSE
#     )
#   }
#   else
#   {
#     CCF_QC = Reduce(dplyr::bind_rows, lapply(x$CCF_estimates, function(x) x$QC_table)) %>%
#       dplyr::select(karyotype, QC) %>%
#       dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
#       dplyr::arrange(karyotype) %>%
#       dplyr::mutate(
#         value = 1,
#         lab.ypos = cumsum(value) - 0.5 * value,
#         QC = paste(QC),
#         label = karyotype,
#         type = 'CCF')
#   }
#
#   QC_table = dplyr::bind_rows(peaks_QC, CCF_QC)
#   QC_table$karyotype = factor(QC_table$karyotype, all_karyptypes)
#   QC_table$type = factor(QC_table$type, levels = c('Peaks', 'CCF'))
#
#   QC_table %>%
#     dplyr::mutate(
#       QC = ifelse(!is.na(QC) & QC == "NA", NA, QC)
#     )
# }
#
# # Page reports pipeline
# page1_pcnaqc = function(x)
# {
#   N = ifelse(x$n_snvs > 3000, 3000, x$n_snvs)
#
#   # Segments
#   A = cowplot::plot_grid(
#     CNAqc::plot_gw_depth(x, N = N) + labs(
#       title = paste0(
#         "Ploidy ",
#         x$ploidy,
#         ', most prevalent karyotype is ',
#         CNAqc:::get_prevalent_karyotype(x),
#         " (Major:minor).\n"
#       )
#     ),
#     CNAqc::plot_segments(x),
#     align = 'v',
#     nrow = 2,
#     rel_heights = c(.35, .8)
#   )
#
#
#   # Histograms plots
#   W = ggpubr::ggarrange(
#     CNAqc::plot_data_histogram(x, which = "VAF"),
#     CNAqc::plot_data_histogram(x, which = "DP"),
#     CNAqc::plot_data_histogram(x, which = "CCF"),
#     ncol = 3,
#     common.legend = T,
#     legend = 'bottom'
#   )
#
#   # Genome distribution and QC
#   E = ggpubr::ggarrange(
#     CNAqc::plot_karyotypes(x),
#     CNAqc::plot_karyotypes(x, type = 'number'),
#     common.legend = TRUE,
#     legend = 'bottom'
#   )
#
#   G = CNAqc::plot_qc(x) + labs(title = paste0("Purity/ ploidy (peaks) and CCF estimates (% assignable)"))
#
#
#   H = cowplot::plot_grid(
#     E,
#     G,
#     ncol = 2,
#     rel_widths = c(1, 1),
#     labels = c("C", 'D'),
#     align = 'h',
#     axis = 'x'
#   )
#
#   ggpubr::ggarrange(
#     A,
#     ggplot() + geom_blank(),
#     W,
#     ggplot() + geom_blank(),
#     H,
#     ncol = 1,
#     heights = c(5.8, 1, 4, 1, 4.5),
#     labels = c("A", '', "B", '', "C")
#   )
# }
#
# page2_pcnaqc = function(x)
# {
#   sc = getOption('CNAqc_cex', default = 1)
#
#   A = cowplot::plot_grid(
#     ggplot() + theme_void(base_size = 10  * sc) + labs(
#       title = paste0(
#         "Purity and ploidy QC via peak detection (purity p = ",
#         x$purity * 100,
#         "%)\n"
#       )
#     ),
#     CNAqc::plot_peaks_analysis(x),
#     align = 'v',
#     nrow = 2,
#     rel_heights = c(.1, .95)
#   )
#
#   N = x$n_snvs * 0.01
#   D = median(x$cna$length) * .1
#
#   segment_ids = x$cna %>%
#     filter(n > !!N, length > D) %>%
#     pull(segment_id)
#
#   if (length(segment_ids) == 0)
#     B = CNAqc:::eplot()
#   else
#   {
#     B = CNAqc::inspect_segment(x, n = N, l = D) +
#       labs(
#         title = paste0(
#           "Segments with >",
#           N,
#           " mutations (1% of total) and >",
#           D,
#           ' bases (10% of median length)\n'
#         ),
#         subtitle = NULL
#       )
#   }
#
#   C1 = CNAqc::plot_arm_fragmentation(x) +
#     labs(title = paste0("Arm-level segments fragmentation (% relative to arm length)\n"))
#   C2 = CNAqc::plot_segment_size_distribution(x) +
#     labs(title = paste0("Segments size distribution\n"))
#
#   cowplot::plot_grid(
#     A,
#     ggplot() + geom_blank(),
#     B,
#     ggplot() + geom_blank(),
#     ggpubr::ggarrange(C1, C2),
#     ncol = 1,
#     rel_heights = c(4, 1, 4, 1, 4),
#     labels = c("A", '', "B", '', "C"),
#     scale = c(.9, 1, 1, 1, 1),
#     align = 'v',
#     axis = 'y'
#   )
#
# }
#
# page3_pcnaqc = function(x)
# {
#   CNAqc::plot_CCF(x) %>%
#     evoverse:::img_annot_t(t = paste0(
#       "QC PASS (per karyotype) if we can assign >90% of mutation burden\n"
#     ))
# }
#
# # Report one-page style for the Data QC pipeline
# report_onepage_cnaqc_pipeline = function(x)
# {
#   # Page one -- segments | vaf | qc
#   # Page two -- Peak analysis and fragmentation
#   # Page three -- Peak analysis and fragmentation
#   first_page = evoverse:::page1_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#   second_page = evoverse:::page2_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#   third_page = evoverse:::page3_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#
#   # Report assembly
#   p1_tit = 'Data quality check (QC) from matched Copy Number and mutation data.'
#   p2_tit = 'Purity/ ploidy QC (A), data (B) and fragmentation (C).'
#   p3_tit = 'Cancer Cell Fractions estimation and QC.'
#
#   # Special report - 1 page wide collated
#   ggpubr::ggarrange(
#       first_page %>% evoverse:::img_annot_t(p1_tit),
#       second_page %>% evoverse:::img_annot_t(p2_tit),
#       third_page %>% evoverse:::img_annot_t(p3_tit),
#       ncol = 3
#     )
# }
#
#
# # Report multi-page style for the Data QC pipeline
# report_multipage_cnaqc_pipeline = function(x, f)
# {
#   # Page one -- segments | vaf | qc
#   # Page two -- Peak analysis and fragmentation
#   # Page three -- Peak analysis and fragmentation
#   first_page = evoverse:::page1_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#   second_page = evoverse:::page2_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#   third_page = evoverse:::page3_pcnaqc(x$fit) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#
#   # Report assembly
#   # tool <-
#   #   magick::image_read("https://caravagn.github.io/CNAqc/reference/figures/logo.png")
#   sample_logo = NA
#
#   # Frontline
#   cp_e_params = function() {
#     B = brewer.pal('Spectral', n = 5)
#
#     # QC
#     l_PASS = cut(x$QC$f_PASS, seq(0, 100, by = 20)) %>% as.numeric
#     qc_col = brewer.pal('Spectral', n = 5)[l_PASS]
#
#     grid::grid.rect(
#       gp = grid::gpar(fill = alpha(qc_col, .3)),
#       vp = grid::viewport(layout.pos.row = 2,
#                           layout.pos.col = 1),
#       width = 0.8,
#       height = 0.2
#     )
#     #
#     grid::grid.text(
#       paste0("PASS rate ", x$QC$f_PASS),
#       gp = grid::gpar(col = qc_col),
#       vp = grid::viewport(layout.pos.row = 2,
#                           layout.pos.col = 1)
#     )
#   }
#
#   p1_tit = 'Data quality check (QC) from matched Copy Number and mutation data.'
#   p2_tit = 'Purity/ ploidy QC (A), data (B) and fragmentation (C).'
#   p3_tit = 'Cancer Cell Fractions estimation and QC.'
#
#   message("Saving multi-page PDF report")
#
#   start_report(size = 'a4', filename = f) %>%
#     add_cover_page(
#       logo = NA,
#       title = x$description,
#       subtitle = p1_tit,
#       logo_size = 0.15,
#       theme = 'basic',
#       creaton_time = paste0(
#         'evoverse pipeline (caravagn.github.io/evoverse), ',
#         Sys.time()
#       ),
#       need_footer = T,
#       extra_layout_params = cp_e_params
#     ) %>%
#     add_new_page(
#       plot = first_page,
#       title = sample,
#       subtitle = 'CNA Segments (A), mutations (B), genome coverage (C) and QC (D).',
#       need_header = TRUE,
#       header_color = "#f44242",
#       logo = sample_logo,
#       logo_size = .35,
#       theme = "flashy",
#       need_footer = TRUE,
#       footer_text = "Copyright © Giulio Caravagna | 2020 | evoverse (caravagn.github.io/evoverse)",
#       footer_color = "#3d3c3c"
#     ) %>%
#     add_new_page(
#       plot = second_page,
#       title = sample,
#       subtitle = p2_tit,
#       need_header = TRUE,
#       header_color = "#f44242",
#       logo = sample_logo,
#       logo_size = .35,
#       theme = "flashy",
#       need_footer = TRUE,
#       footer_text = "Copyright © Giulio Caravagna | 2020 | evoverse (caravagn.github.io/evoverse)",
#       footer_color = "#3d3c3c"
#     ) %>%
#     add_new_page(
#       plot = third_page,
#       title = sample,
#       subtitle = p3_tit,
#       need_header = TRUE,
#       header_color = "#f44242",
#       logo = sample_logo,
#       logo_size = .35,
#       theme = "flashy",
#       need_footer = TRUE,
#       footer_text = "evoverse GitHub Pages: caravagn.github.io/evoverse | @gcaravagna | 2020",
#       footer_color = "#3d3c3c"
#     ) %>%
#     end_report()
# }
