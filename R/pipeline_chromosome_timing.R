#
# cna = CNAqc::example_dataset_CNAqc$cna
# mutations = CNAqc::example_dataset_CNAqc$snvs
# purity = CNAqc::example_dataset_CNAqc$purity

# x = chromosome_timing_pipeline(mutations, cna = cna, purity = purity, auto_setup = 'FAST')

#' Pipeline to time aneuploidy with MOBSTER.
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
#' x = pipeline_chromosome_timing(mutations, cna = cna, purity = purity, auto_setup = 'FAST')
#' print(x)
pipeline_chromosome_timing = function(mutations,
                                      cna = NULL,
                                      purity = NULL,
                                      timeable = c('2:0', '2:1', '2:2'),
                                      min_muts = 50,
                                      ...
                                      )
{
  pio::pioHdr("Evoverse", italic('Copy Number timing pipeline'))
  cat('\n')

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
    cli::cli_process_start("Using CNA data to subset mutations.")
    cat("\n")

    cna_obj = CNAqc::init(mutations, cna, purity)
    mutations = cna_obj$snvs

    available_karyo = cna_obj$n_karyotype[timeable]
    available_karyo = available_karyo[!is.na(available_karyo)]
    timeable = names(available_karyo)

    cat("\n")
    cli::cli_process_done()
  }

  # Deconvolution function
  mfits = deconvolution_mobster_karyotypes(
    mutations = mutations,
    karyotypes = timeable,
    min_muts = min_muts,
    ...
  )

  # Assemble tables, plots and perform QC
  results = wrap_up_pipeline_mobster(mfits, qc_type = "R", cna_obj)
  results$mobster = mfits
  results$input = list(mutations = mutations, cna = cna, purity = purity)


  return(results)
}

