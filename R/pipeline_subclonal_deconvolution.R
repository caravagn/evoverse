#
# cna = CNAqc::example_dataset_CNAqc$cna
# mutations = CNAqc::example_dataset_CNAqc$snvs
# purity = CNAqc::example_dataset_CNAqc$purity

# target_k = c("1:0", '1:1', '2:1', '2:2')
# target_k = target_k[cna_obj$ploidy]
# x = chromosome_timing_pipeline(mutations, cna = cna, purity = purity, auto_setup = 'FAST')

#' Pipeline to perform subclonal deconvolution with MOBSTER.
#'
#' @description
#'
#' For a list of input mutations, CNA and purity data from a tumour sample,
#' MOBSTER fits are computed to determine proportions of mutations in each
#' aneuploidy state (copy state). Assignments and plots are returned, and
#' an automatic QC pipeline (supervised logistic classifier) is used to
#' flag which fits pass quality check. The trained set for the classifier
#' consists in n = 377 colorectal whole genome samples.
#'
#' @param mutations Input mutations in the format for package \code{CNAqc}.
#' @param cna Absolute Copy Number Alterations in the format for package \code{CNAqc}.
#' @param purity Tumour purity, in [0, 1].
#' @param timeable Karyotypes that can be timed, expressed as "Major:minor" notation. Default
#' are the following states \code{timeable = c('2:0', '2:1', '2:2')}.
#' @param min_muts Skip analysing karyotypes with less than \code{min_muts} mutations.
#' @param ... Parameters passed \code{mobster_fit} in package \code{mobster}.
#'
#' @return
#' @export
#'
#' @examples
#' # We use data released with the CNAqc package
#'
#' cna = CNAqc::example_dataset_CNAqc$cna
#' mutations = CNAqc::example_dataset_CNAqc$snvs
#' purity = CNAqc::example_dataset_CNAqc$purity
#'
#' x = pipeline_subclonal_deconvolution(mutations, cna = cna, purity = purity, auto_setup = 'FAST')
#' print(x)
pipeline_subclonal_deconvolution = function(mutations,
                                      cna = NULL,
                                      purity = NULL,
                                      karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                      CCF_karyotypes = karyotypes,
                                      min_muts = 50,
                                      description = "MOBSTER deconvolution",
                                      ...
)
{

  pio::pioHdr("Evoverse", italic('Subclonal deconvolution pipeline'))
  cat('\n')

  #
  # 1) Check input data, subset by CNA data if available
  #

  # x must be a dataset for MOBSTER
  if(!is.data.frame(mutations) & !is.matrix(mutations))
    stop("mutations: Must use a dataset for MOBSTER!")

  if(!is.null(cna) & !is.data.frame(cna) & !is.matrix(cna))
    stop("cna: Must use a dataset for MOBSTER!")

  if(!is.null(purity) & (purity > 1 | purity <= 0))
    stop("Purity must be in 0/1.")

  # Apply CNA mapping and retain only mappable mutations
  cna_obj = NULL
  if(!is.null(cna))
  {
    cli::cli_h1("Using CNA data to subset mutations.")
    cat("\n")

    cna_obj = CNAqc::init(mutations, cna, purity)
    mutations = cna_obj$snvs
  }

  #
  # 2) MOBSTER analysis of karyotypes and CCF if required
  #
  mfits = deconvolution_mobster_karyotypes(
      mutations = mutations,
      karyotypes = karyotypes,
      min_muts = min_muts,
      ...
    )

  if(!is.null(CCF_karyotypes) | !is.na(CCF_karyotypes))
  {
    if(is.null(cna_obj)) stop("You asked for CCF analysis but did not provide CNA data (cna = NULL), cannot compute.")

    CCF_fit = deconvolution_mobster_CCF(cna_obj, CCF_karyotypes = CCF_karyotypes, min_muts = min_muts)
    cna_obj = CCF_fit$cna_obj

    mfits = append(mfits, list(`CCF` = CCF_fit$fits))
  }

  #
  # 3) Wrap MOBSTER results for all this fits, apply QC
  #
  results =  wrap_up_pipeline_mobster(mfits,
                                     qc_type = "D",
                                     cna_obj,
                                     karyotypes = karyotypes,
                                     add_CCF = ("CCF" %in% names(mfits))
                                     )
  results$mobster = mfits
  results$input = list(mutations = mutations, cna = cna, purity = purity, cnaqc_obj = cna_obj)

  #
  # 4) BMix -- Binomial model -- on each karyotype non-tail mutations
  #
  runner = function(k)
  {
    if(all(is.null(results$mobster[[k]]))) return(NULL)

    non_tail = mobster::Clusters(results$mobster[[k]]$best) %>%
      dplyr::filter(cluster != "Tail") %>%
      dplyr::select(NV, DP) %>%
      data.frame

    Kbeta = min(results$mobster[[k]]$best$Kbeta * 2, 4)

    fit_readcounts = BMix::bmixfit(
      non_tail,
      K.BetaBinomials = 0,
      K.Binomials = 1:Kbeta,
      samples = 3)

    fit_readcounts$input = non_tail
    fit_readcounts
  }

  bmix_fits = lapply(names(results$mobster), runner)
  names(bmix_fits) = names(results$mobster)

  results$bmix = bmix_fits

  # Set a panel like the one above, same dimension
  bmix_panel = lapply(bmix_fits,
                      function(x)
                      {
                        if (all(is.null(x)))
                          return(ggplot() + geom_blank())
                        BMix::plot_clusters(x, x$input)
                      })

  bmix_panel = ggarrange(plotlist = append(list(ggplot() + geom_blank()), bmix_panel),
                         nrow = 1,
                         ncol = length(bmix_panel) + 1,
                         labels = names(results$mobster))

  #
  # 5) Assemble a one-page final figure with both MOBSTER and BMix panels
  #
  figure =  ggpubr::ggarrange(
    results$figure,
    bmix_panel,
    nrow = 2,
    ncol = 1
    )

  caption = paste0("", Sys.time(), '. evoverse pipeline for subclonal deconvolution. QC: ', results$qc$QC_type[1])

  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(bquote(bold("Dataset. ") ~ description), hjust = 0, x = 0),
    bottom = ggpubr::text_grob(bquote(caption, hjust = 0, x = 0, size = 8))

  results$figure = figure
  results$description = description

  #
  # 6) Generate tables
  #

  # QC table for mobster
  results$tables$qc_mobster = results$qc
  results$qc = NULL

  # Clustering assignments table contain both MOBSTER and BMix results for the same sample
  results$tables$mobster_assignments = results$assignments

  bmix_assignments = lapply(
    names(bmix_fits),
    function(x)
    {
      if(all(is.null(bmix_fits[[x]]))) return(NULL)
      x = BMix::to_string(bmix_fits[[x]])
      colnames(x) = paste0("BMix_", colnames(x))
      x
    })


  results$assignments = results$qc = NULL

  return(results)
}

