#' Plot a dataset (complex figure).
#'
#' @description Create a complex figure via \code{ggpubr} assembling several plots. The figure
#' reports 1) the adjusted VAF (full spectrum), 2) the zoomed version of the adjusted
#' VAF for values below 0.20 (low-frequency spectrum), 3) the depth against the VAF
#' of each mutation and 4) the depth distribution.
#'
#'
#' @param x A dataset of class \code{"mbst_data"}.
#' @param ...
#'
#' @return A ggpubr figure object.
#' @export
#'
#' @examples
#' TODO
plot.mbst_data = function(x, ...)
{
  check_is_mobster_mvdata(x)

  params = as.list(substitute(list(...)))[-1L]
  which_cluster = params['clusters']

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Plots without clusters
  # =-=-=-=-=-=-=-=-=-=-=-=-
  if(which_cluster != 'MOBSTER' & which_cluster != 'VIBER')
  {
    message("Printing data - use cluster = 'MOBSTER' or cluster = 'VIBER' to plot clusters")

    # 2D case, without clusters
    if (length(x$samples) > 1)
    {
      pairs = combn(x$samples, 2)

      figures = apply(pairs, 2, function(w)
      {
        plot_2D_VAF(x, w[1], w[2], ...)
      })
    }
    else
    {
      figures = plot_1D_VAF(x, x$samples[1])
    }
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Plots with MOBSTER cluster
  # =-=-=-=-=-=-=-=-=-=-=-=-
  if(which_cluster == 'MOBSTER')
  {
    message("Printing data augmented with MOBSTER clusters")

    # 2D case, without clusters
    if (length(x$samples) > 1)
    {
      pairs = combn(x$samples, 2)

      figures = apply(pairs, 2, function(w)
      {
        plot_2D_VAF_MOBSTER(x, w[1], w[2], ...)
      })
    }
    else
    {
      figures = plot(x$fit_MOBSTER[[1]]$best)
    }

  }

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Plots with VIBER cluster
  # =-=-=-=-=-=-=-=-=-=-=-=-
  if(which_cluster == 'VIBER')
  {
    message("Printing data augmented with VIBER clusters - using plotting functions from the VIBER package.")

    figures = plot(x$fit_VIBER)
  }


  return(figures)
}

#
