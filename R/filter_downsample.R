#' Downample data.
#'
#' @description Randomly downsample the data contained in a dataset to a
#' certain number of mutations.
#'
#' @param x A dataset of class \code{"mbst_data"}.
#' @param N The number of mutations to retain. If \code{x} contains less
#' than \code{N} mutations, \code{x} is returned.
#'
#' @return A dataset with at most \code{N} mutations.
#' @export
#'
#' @examples
#' TODO
filter_downsample = function(x, n)
{
  check_is_mobster_mvdata(x)
  stopifnot(n > 0)

  if(n < N(x))
  {
    pioStr(paste0("Subsampling ", n, "entries out of ", N(x)))

    k = keys(x)
    ids = sample(k, length(k) - n)

    x = delete_entries(x, ids)
    x = logOp(x, paste0("Subsampled to N = ", n, ""))
  }

  x
}
