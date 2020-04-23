#' Return the number of mutations in the dataset.
#'
#' @description Return the number of mutations in the dataset.
#'
#' @param x An `evoverse` object.
#'
#' @return The number of mutations in the dataset.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' N(example_evoverse)
N = function(x) {
  x$mutations_locations %>% nrow
}

#' Return the number of samples in the dataset.
#'
#' @description Return the number of samples in the dataset.
#'
#' @param x An `evoverse` object.
#'
#' @return The number of samples in the dataset.
#' @export
#'
#' @examples
#' data('example_evoverse')
#' S(example_evoverse)
S = function(x) {
  length(x$samples)
}


