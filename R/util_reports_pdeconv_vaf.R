# # First page of the report for the pipeline for subclonal deconvolution from raw VAF data
# page1_pdeconv_vaf = function(x)
# {
#   mobster_fits = x$mobster
#   bmix_fits = x$bmix
#   cna_obj = x$input$cnaqc
#   qc_table = x$table$summary
#
#
#   USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
#
#   groups = names(x$mobster)
#   empty_panel = CNAqc:::eplot()
#
#   A = B = C = D = E = empty_panel
#
#   # A genome data ~ CNA plot
#   if (!is.null(cna_obj))
#     # Segments
#     A = cowplot::plot_grid(
#       CNAqc::plot_depth(cna_obj, N = 3000) + labs(
#         title = paste0(
#           "Ploidy ",
#           cna_obj$ploidy,
#           ', most prevalent karyotype is ',
#           CNAqc:::get_prevalent_karyotype(cna_obj),
#           " (Major:minor).\n"
#         )
#       ),
#       CNAqc::plot_segments(cna_obj),
#       align = 'v',
#       nrow = 2,
#       rel_heights = c(.35, .8)
#     )
#
#   # B - Karyotype percentage
#   B = CNAqc::plot_karyotypes(cna_obj) + labs(title = "Karyotypes coverage\n")
#
#   # C button plot
#   expanded = x$table$QC_table %>%
#     dplyr::rename(deconvolution = QC, peaks = QC_peaks) %>%
#     dplyr::select(karyotype, peaks, deconvolution) %>%
#     reshape2::melt(id = 'karyotype')
#
#   Nt = sum(x$table$QC_table$N, na.rm = T)
#
#   expanded = expanded %>%
#     left_join(x$table$QC_table %>%
#                 dplyr::select(-starts_with("QC")) %>%
#                 dplyr::mutate(N = log(N / Nt)),
#               by = 'karyotype') %>%
#     mutate(N = ifelse(is.na(N), 1, N))
#
#   C = expanded %>%
#     ggplot() +
#     geom_tile(
#       aes(
#         x = karyotype,
#         y = variable,
#         width = .8,
#         height = .8,
#         color = paste(architecture),
#         # size = as.numeric(N),
#         fill = paste(value)
#       ),
#       size = 2
#     ) +
#     labs(title = "Summary QC results\n",
#          y = 'QC result (NA, Not available)') +
#     scale_fill_manual(values = c(
#       `PASS` = 'forestgreen',
#       `FAIL` = 'indianred3',
#       `NA` = 'gainsboro'
#     )) +
#     CNAqc:::my_ggplot_theme() +
#     guides(fill = guide_legend('QC'),
#            color = guide_legend('Architecture', override.aes = list(fill = NA))) +
#     scale_color_manual(values = c(
#       `NA` = NA,
#       `Monoclonal` = 'orange',
#       `Polyclonal` = 'purple4'
#     ))
#
#   # B = gridExtra::tableGrob(x$table$QC_table)
#
#
#
#
#   # Assemble a one-page final figure with both MOBSTER and BMix panels
#   ggpubr::ggarrange(
#     A,
#     ggplot() + geom_blank(),
#     cowplot::plot_grid(
#       B,
#       C,
#       rel_widths = c(1, 2),
#       nrow = 1,
#       align = 'h',
#       labels = c("B", "C")
#     ),
#     ncol = 1,
#     heights = c(4, 1, 4),
#     labels = c("A", '', '', "")
#   )
#
#   # # Set the figure title and captions
#   # figure = ggpubr::annotate_figure(
#   #   figure,
#   #   top = ggpubr::text_grob(bquote(bold("Dataset. ") ~ .(x$description)), hjust = 0, x = 0, size = 15),
#   #   bottom = ggpubr::text_grob(bquote(.(x$log)), hjust = 0, x = 0, size = 10)
#   # )
#
#
# }
#
# # Second page of the report for the pipeline for subclonal deconvolution from raw VAF data
# page2_pdeconv_vaf = function(x)
# {
#   mobster_fits = x$mobster
#   bmix_fits = x$bmix
#   cna_obj = x$input$cnaqc
#   qc_table = x$table$summary
#
#   USE_KARYOTYPES = c("1:0", '1:1', '2:0', '2:1', '2:2')
#
#   groups = names(mobster_fits)
#   empty_panel = CNAqc:::eplot()
#
#   D = E = empty_panel
#
#   # MOBSTER plots, sourrounded by a coloured box by QC
#   D = lapply(groups,
#              function(y)
#              {
#                if (all(is.null(mobster_fits[[y]])))
#                  return(CNAqc:::eplot())
#
#                qc_entries = qc_table %>% dplyr::filter(karyotype == !!y)
#
#                qc = ifelse(qc_entries$QC == "FAIL", "indianred3", 'forestgreen')
#
#                w = mobster_fits[[y]]
#                subl = paste0('n = ',
#                              w$N,
#                              ", ",
#                              paste(
#                                names(w$pi),
#                                ' ',
#                                round(w$pi * 100),
#                                '%',
#                                collapse = ', ',
#                                sep = ''
#                              ))
#
#                mobster::plot.dbpmm(mobster_fits[[y]]) +
#                  labs(
#                    caption = NULL,
#                    subtitle = NULL,
#                    title = bquote(.(y) ~ " (" * p["PASS"] ~ '=' ~ .(qc_entries$QC_prob) *
#                                     ');' ~ .(subl))
#                  ) +
#                  theme(title = element_text(color = qc),
#                        panel.border = element_rect(colour = qc,
#                                                    fill = NA))
#              })
#
#   # If there is a second panel, we put it below
#   if (!all(is.null(bmix_fits)))
#   {
#     # BMIx: a panel like the one above, same dimension
#     E = lapply(groups,
#                function(x)
#                {
#                  if (all(is.null(bmix_fits[[x]])))
#                    return(empty_panel)
#
#                  w = bmix_fits[[x]]
#
#                  subl = paste0('n = ',
#                                nrow(w$input),
#                                ", ",
#                                paste(
#                                  names(w$pi),
#                                  ' ',
#                                  round(w$pi * 100),
#                                  '%',
#                                  collapse = ', ',
#                                  sep = ''
#                                ))
#
#
#                  BMix::plot_clusters(bmix_fits[[x]], bmix_fits[[x]]$input %>% dplyr::select(NV, DP)) +
#                    scale_fill_brewer(palette = 'Set2') +
#                    labs(title = subl,
#                         subtitle = NULL)
#                })
#   }
#
#   inactive_plots = sapply(mobster_fits, function(w)
#     all(is.null(w)))
#   squeeze = rep(1, length(inactive_plots))
#   squeeze[!inactive_plots] = 2
#
#   D = ggpubr::ggarrange(
#     plotlist = D,
#     ncol = 1,
#     nrow = length(groups),
#     heights = squeeze
#   )
#
#   E = ggarrange(
#     plotlist = E,
#     ncol = 1,
#     nrow = length(groups),
#     heights = squeeze
#   )
#
#   DE = ggpubr::ggarrange(D, E, ncol = 2)
#
#   return(DE %>% evoverse:::img_annot_t("mobster and BMix fits per karyotype"))
# }
#
# # Report multi-page style for the Data QC pipeline
# report_multipage_deconv_pipeline = function(x, f, sample, collate = F)
# {
#   require(myprettyreport)
#
#   # QC
#   PASS = x$table$QC == "PASS"
#   qc_col = ifelse(PASS, 'forestgreen', 'indianred3')
#
#   # Page one -- CNA | genome | QC
#   first_page = page1_pdeconv_vaf(x) %>% evoverse:::img_add_margin(fact = 1, w = 1, h = 2)
#   second_page = page2_pdeconv_vaf(x) %>% evoverse:::img_add_margin(fact = .5, w = 1, h = 1)
#
#   # tool <-
#   #   magick::image_read("https://caravagn.github.io/CNAqc/reference/figures/logo.png")
#   tool = NA
#   sample_logo = tool
#
#   # Frontline
#   cp_e_params = function() {
#     grid::grid.rect(
#       gp = grid::gpar(fill = alpha(qc_col, .3)),
#       vp = grid::viewport(layout.pos.row = 2,
#                           layout.pos.col = 1),
#       width = 0.8,
#       height = 0.2
#     )
#     #
#     col_text = ifelse(x$table$architecture == "Monoclonal", 'orange', 'purple4')
#     if (is.na(col_text))
#       col_text = 'black'
#
#     grid::grid.text(
#       paste0(x$table$architecture, ' architecture ~ QC ', x$table$QC, '.'),
#       gp = grid::gpar(col = col_text),
#       vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1)
#     )
#   }
#
#   p1_tit = 'CNA segments (A), karyotype coverage (B) and summary QC (C).'
#   p2_tit = 'Subclonal deconvolution from raw VAF data, per karyotype.'
#
#   # Special report - 1 page wide collated
#   if (collate == TRUE)
#   {
#     message("Saving single-page collated PDF report")
#
#     ggpubr::ggarrange(
#       first_page %>% evoverse:::img_annot_t(p1_tit),
#       second_page %>% evoverse:::img_annot_t(p2_tit),
#     ) %>%
#       ggsave(
#         filename = f,
#         height = 250,
#         width = 0.7 * 297,
#         units = 'mm'
#       )
#   }
#   else
#   {
#     message("Saving multi-page PDF report")
#
#     start_report(size = 'a4', filename = f) %>%
#       add_cover_page(
#         logo = sample_logo,
#         title = sample,
#         subtitle = 'Evoverse pipeline for subclonal deconvolution from raw VAF data, by karyotype.',
#         logo_size = 0.15,
#         theme = 'basic',
#         creaton_time = paste0(
#           'evoverse pipeline (caravagn.github.io/evoverse), ',
#           Sys.time()
#         ),
#         need_footer = T,
#         extra_layout_params = cp_e_params
#       ) %>%
#       add_new_page(
#         plot = first_page,
#         title = sample,
#         subtitle = p1_tit,
#         need_header = TRUE,
#         header_color = "#f44242",
#         logo = sample_logo,
#         logo_size = .35,
#         theme = "flashy",
#         need_footer = TRUE,
#         footer_text = "evoverse GitHub Pages: caravagn.github.io/evoverse | @gcaravagna | 2020",
#         footer_color = "#3d3c3c"
#       ) %>%
#       add_new_page(
#         plot = second_page,
#         title = sample,
#         subtitle = p2_tit,
#         need_header = TRUE,
#         header_color = "#f44242",
#         logo = sample_logo,
#         logo_size = .35,
#         theme = "flashy",
#         need_footer = TRUE,
#         footer_text = "evoverse GitHub Pages: caravagn.github.io/evoverse | @gcaravagna | 2020",
#         footer_color = "#3d3c3c"
#       ) %>%
#       end_report()
#   }
#
# }
