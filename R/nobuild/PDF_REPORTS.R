#' Title
#'
#' @param x
#' @param file
#' @param cex
#'
#' @return
#'
#' @import myprettyreport
#' @import patchwork
#'
#' @export
#'
#' @examples
pdf_report = function(x, file, cex = 1)
{
  require(patchwork)
  require(myprettyreport)

  # Set global cex value for all plots
  options(CNAqc_cex = cex)

  # Running
  cli::cli_h1("PDF report assembly")
  cli::cli_process_start("Creating report in file {.field {file}}")
  cat('\n')

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wrapper function for evopipe_qc
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(inherits(x, "evopipe_qc"))
  {
    options(CNAqc_cex = cex)

    # Figure assembly
    suppressWarnings(
      evoverse:::report_multipage_cnaqc_pipeline(x, f = file)
    )

  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wrapper function for evopipe_qc
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cli::cli_process_done()
}
