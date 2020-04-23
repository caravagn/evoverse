#' Randomly downsample a dataset.
#'
#' @description Randomly downsample a dataset to contain a certain number of mutations.
#' If the dataset contains less than the requirednumber of mutations, the input object is returned.
#'
#' @param x An `evoverse` object.
#' @param n The number of mutations to retain. If \code{x} contains less
#' than \code{n} mutations, \code{x} is returned.
#'
#' @return A dataset with at most \code{n} mutations.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' N(example_evoverse)
#'
#' N(filter_downsample(example_evoverse, n = 100))
filter_downsample = function(x, n)
{
  check_is_mobster_mvdata(x)
  stopifnot(n > 0)

  if(n < N(x))
  {
    pioStr("Subsampling mutations:", n, "out of", N(x), suffix = '\n')

    k = keys(x)
    ids = sample(k, n)

    x = delete_entries(x, ids = setdiff(k, ids))
    x = logOp(x, paste0("Subsampled to N = ", n, ""))
  }

  x
}
