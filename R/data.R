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

#' Example data in mvMOBSTER format.
#'
#' @docType data
#'
#' @usage data(example_mvmobster)
#'
#' @format A tibble that represents the cytobands for hg19 chromosome
#' (chr, from, to, cytoband and region).
#'
#' @keywords datasets
#'
#' @references hg19
#'
#' @examples
#' data(example_mvmobster)
#' head(example_mvmobster)
#'
"example_mvmobster"

#' Example data in the wide and long formts supported by mvMOBSTER
#'
#' @docType data
#'
#' @usage data(example_input_formats)
#'
#' @format 4 tibbles, wide/long and mutations/copy number that
#' show the format required to create a mvMOBSTER object with
#' function \link{\code{dataset}}.
#'
#' @keywords datasets
#'#'
#' @examples
#' data(example_input_formats)
#' head(example_input_formats)
"example_input_formats"
