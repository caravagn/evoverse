#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.evopipe_deconv = function(x)
{
  stopifnot(inherits(x, 'evopipe_deconv'))

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ Evoverse ] {.value {x$description}}")),
      '{.value {x$type}}'
    )
  )

  # Print pipeline objects
  # cli::cli
  # cat(cli::rule(line = '--&<--'))
  cat("\n")
  print(x$input$cnaqc)

  cat("\n")

  cli::cli_alert_info("This tumour is {.field {x$table$architecture}}.")
  cat("\n")

  df = x$table$summary
  df = df[, names(df) %in% c('karyotype', 'N', 'K_beta', 'tail', 'Shape_Tail', 'QC', 'QC_prob', 'BMix_K')]

  print(x$table$QC_table)
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.evopipe_deconv = function(x)
{
  stopifnot(inherits(x, 'evopipe_deconv'))

  mobster_fits = x$mobster
  bmix_fits = x$bmix
  cna_obj = x$input$cnaqc

  qc_table = x$table$summary

  # 2 x 1 CNAqc plot
  groups = names(mobster_fits)
  empty_panel = CNAqc:::eplot()

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(
    groups,
    function(y)
    {
      if (all(is.null(mobster_fits[[y]])))
        return(CNAqc:::eplot())

      qc_entries = qc_table %>% dplyr::filter(karyotype == !!y)

      qc = ifelse(qc_entries$QC == "FAIL", "indianred3", 'forestgreen')

      mobster::plot.dbpmm(mobster_fits[[y]]) +
        labs(title = bquote("QC " ~ .(qc_entries$QC) ~ p["PASS"] ~
                              '=' ~ .(qc_entries$QC_prob))) +
        theme(title = element_text(color = qc),
              panel.border = element_rect(
                colour = qc,
                fill = NA,
                size = 5
              ))
  })


  # CNA plot
  cna_plot = empty_panel
  if(!is.null(cna_obj)) cna_plot = CNAqc::plot_segments(cna_obj, circular = TRUE)

  # Top panel: CNA + MOBSTER
  mob_fits_plot =  append(list(cna_plot), mob_fits_plot)
  figure = ggpubr::ggarrange(
    plotlist = mob_fits_plot,
    nrow = 1,
    ncol = length(groups) + 1,
    labels = c("CNA", groups)
  )

  # If there is a second panel, we put it below
  if(!all(is.null(bmix_fits)))
  {
    # BMIx: a panel like the one above, same dimension
    bmix_panel = lapply(
      groups,
      function(x)
      {
        if (all(is.null(bmix_fits[[x]]))) return(empty_panel)
        BMix::plot_clusters(bmix_fits[[x]], bmix_fits[[x]]$input %>% dplyr::select(NV, DP))
      })

    bmix_panel = ggarrange(plotlist = append(list(empty_panel), bmix_panel),
                           nrow = 1,
                           ncol = length(groups) + 1,
                           labels = c("", groups))

    # Assemble a one-page final figure with both MOBSTER and BMix panels
    figure =  ggpubr::ggarrange(
      figure,
      bmix_panel,
      nrow = 2,
      ncol = 1
    )
  }

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(bquote(bold("Dataset. ") ~ .(x$description)), hjust = 0, x = 0, size = 15),
    bottom = ggpubr::text_grob(bquote(.(x$log)), hjust = 0, x = 0, size = 10)
  )

  return(figure)

}
