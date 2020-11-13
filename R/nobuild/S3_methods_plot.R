#' Plot an \code{evoverse} dataset.
#'
#' @description
#'
#' This function can return different types of lists of plots for an \code{evoverse} dataset;
#' the control on which type of plot is generated is given by the ellipsis parameters.
#'
#' In particular, this function returns:
#'
#' \itemize{
#' \item 2D VAF plots (sample versus sample), if there are more than one sample in this dataset.
#' \item or a single histogram plot when there one sample in the dataset.
#' }
#'
#' The plotting style depends on the ellipsis parameters:
#'
#' \itemize{
#' \item without parameters, the mutations are left black. These plots are generated
#' by ether function \link{plot_2D_VAF}, or \link{plot_1D_VAF};
#' \item with \code{clusters = "MOBSTER"}, the mutations are coloured according to the
#' \code{MOBSTER} clusters. These plots in 2D are generated
#' by function \link{plot_2D_VAF_MOBSTER}, and in 1D by the
#' \code{MOBSTER} S3 \code{plot} function (see MOBSTER's documentation);
#' \item with \code{clusters = "VIBER"}, the mutations are plot according to the
#' \code{VIBER} S3 \code{plot} function (see VIBER's documentation);
#' }
#'
#' @param x An  \code{evoverse} object.
#' @param ... Only one extra parameter is supported: \code{clusters}, which can be
#' set to \code{"MOBSTER"} or \code{"VIBER"}.
#'
#' @return A list of \code{ggplot} figures.
#'
#' @seealso Plotting functions used by this S3 function are \link{\code{plot_2D_VAF}},
#' \link{\code{plot_2D_VAF_MOBSTER}} and \link{\code{plot_1D_VAF}}.
#'
#' @export
#'
#' @examples
#' data('example_evoverse')
#'
#' # Plot the data
#' plot(example_evoverse)
#'
#' # Plot the data augmented with MOBSTER clusters
#' plot(example_evoverse, clusters = "MOBSTER")
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
    message("Plotting plain data - use clusters = {'MOBSTER', 'VIBER'} to plot clusters")

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
    message("Plotting data augmented with MOBSTER clusters")

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
    message("Plotting VIBER clusters- using plot function from the VIBER package.")

    figures = plot(x$fit_VIBER)
  }


  return(figures)
}

#
