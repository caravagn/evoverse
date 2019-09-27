#' Coordinates for hg19 chromosomes.
#'
#' @docType data
#'
#' @usage data(chr_coordinates_hg19)
#'
#' @format A tibble that represents the coorinates for hg19 chromosome
#' (chr, length, from, to, centromerStart and centromerEnd)
#'
#' @keywords datasets
#'
#' @references hg19
#'
#' @examples
#' data(chr_coordinates_hg19)
#' chr_coordinates_hg19
#'
"chr_coordinates_hg19"

#' Example \code{evoverse} dataset.
#'
#' @docType data
#'
#' @usage data(example_evoverse)
#'
#' @format A simple dataset created and analyzed with the
#' \code{evoverse} package.
#'
#' @keywords datasets
#'
#' @examples
#' data(example_evoverse)
#' head(example_evoverse)
#'
"example_evoverse"

#' Example data in the wide and long formts supported by the
#' \code{evoverse} package.
#'
#' @docType data
#'
#' @usage data(example_input_formats)
#'
#' @format 4 tibbles, wide/long and mutations/copy number that
#' show the format required to create an  \code{evoverse} object with
#' function \link{\code{dataset}}.
#'
#' @keywords datasets
#'
#' @examples
#' data(example_input_formats)
#' head(example_input_formats)
"example_input_formats"
