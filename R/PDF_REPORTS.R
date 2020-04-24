#' Title
#'
#' @param x
#' @param file
#' @param cex
#'
#' @return
#' @export
#'
#' @examples
pdf_report = function(x, file, cex = 1)
{
  cli::cli_h1("PDF report assembly")
  cli::cli_process_start("Creating report in file {.field {output}}")
  cat('\n')

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wrapper function for evopipe_qc
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(inherits(x, "evopipe_qc"))
  {
    options(CNAqc_cex = cex)

    # Figure assembly
    suppressWarnings(
      evoverse:::report_multipage_cnaqc_pipeline(
        x$fit,
        f = file,
        cex = cex,
        sample = x$description,
        collate = FALSE,
        score = x$QC$f_PASS
      )
    )

    }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wrapper function for evopipe_qc
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cli::cli_process_done()
}
