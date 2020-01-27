.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c('tidyverse', 'pio', 'crayon', 'CNAqc', 'mobster', 'ctree', 'mtree', 'revolver')

  suppressMessages(
    suppressWarnings(
      sapply(requirements, require, character.only = TRUE)
    )
  )

  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)

  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-

  evoverse_welcome_message =  getOption('evoverse_welcome_message', default = TRUE)

  if(evoverse_welcome_message)
  {
    pk = 'evoverse'
    pk_l = 'Cancer Evolution analysis from bulk DNA sequencinng'
    www = "https://caravagn.github.io/evoverse/"
    em = "gcaravagn@gmail.com"

    cli::cli_alert_success(
      'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )

    options(evoverse_welcome_message = FALSE)
  }
}
