#' Extract clustering results.
#'
#' @description This function extracts the data included into the dataset object by using
#' function \code{Data_table}, and then attempts extracting clustering from both MOBSTER
#' and VIBER analyses. If none of these are found, a warning is raised and data are returned.
#' Otherwise, the data tibble is augmented to have one column for each one of the available
#' clusters. If a sample is named `X`, column `X.MOBSTER_clusters` reports the
#' clustering assigments from a MOBSTER analysis of the sample. Instead, VIBER clusters
#' (multivariate), are reported as a unique column `VIBER_clusters`.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#'
#' @return Different types of clustering assigments for every mutation.
#'
#' @export
#'
#' @examples
#' TODO
Clusters = function(x)
{
  check_is_mobster_mvdata(x)

  dtable = Data_table(x)

  if (!has_mobster_fits(x) & !has_viber_fits(x))
  {
    warning("There are no clusters avaialble for this dataset, will return just the data.")
  }

  # MOBSTER clusters extraction
  MOBSTER_clusters = NULL

  if (has_mobster_fits(x))
  {
    list.best = lapply(x$fit_MOBSTER,
                       function(w)
                         return(w$best$data %>% select(-sample, -VAF,-karyotype)))

    MOBSTER_clusters = list.best %>% purrr::reduce(full_join, by = "id")
    colnames(MOBSTER_clusters)[2:ncol(MOBSTER_clusters)] = paste0(names(x$fit_MOBSTER), '.MOBSTER_cluster')

    dtable = dtable %>%
      full_join(MOBSTER_clusters, by = 'id')
  }

  # VIBER clusters extraction
  VIBER_clusters = NULL

  if (has_viber_fits(x))
  {
    VIBER_clusters = bind_cols(x$fit_VIBER$labels, x$fit_VIBER$fit_data)

    dtable = dtable %>%
      full_join(VIBER_clusters, by = 'id')
  }

  dtable
}


