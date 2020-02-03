print.evopipe_deconv = function(x)
{
  stopifnot(inherits(x, 'evopipe_deconv'))

  mobster:::print.dbpmm

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

  x$table$summary %>%
    dplyr::select(karyotype, N, K_beta, tail, Shape_Tail, QC, QC_prob, BMix_K) %>%
    print



}
