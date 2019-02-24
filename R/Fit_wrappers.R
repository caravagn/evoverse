#' Compute MOBSTER fits for the available samples.
#'
#' @description
#'
#' This function wraps a call to the \code{mobster_fit} function
#' from package \code{mobster}, which implements the MOBSTER
#' statistical model. This is a combination of one Pareto power
#' law for a tumour tail, plus a mixture of Beta distributions.
#' The function applies the fitting function to the samples stored in the
#' dataset \code{x}, one by one, and the returned object contains the
#' MOBSTER fits which can be inspected with the functions available
#' in the \code{mobster} package.
#'
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param samples The samples on which MOBSTER will be run.
#' @param ... Extra parameters that will be forwarded to
#' \code{mobster_fit} function from package \code{mobster}.
#'
#' @return A new object with MOBSTER fits for the required
#' samples in \code{x$fit.MOBSTER}.
#'
#' @import mobster
#'
#' @export
#'
#' @examples
#' TODO
mobster_fit_multivariate = function(x, samples = x$samples, ...)
{
  require(mobster)

  x$fit.MOBSTER = lapply(samples,
                         function(w)
                         {
                           pio::pioTit(paste("Fitting with MOBSTER sample", w))

                           data = VAF(x, samples = w) %>%
                             filter(value > 0) %>%
                             spread(variable, value)

                           print(data)

                           mf = mobster::mobster_fit(
                             x = data,
                             ...)

                           mf
                         })

  names(x$fit.MOBSTER) = samples

  cat("\n")
  pio::pioTit("MOBSTER fits")
  for(s in samples) print(x$fit.MOBSTER[[s]]$best)

  # Log update
  x = logOp(x, paste0("Fit MOBSTER to ", paste0(samples, collapse = ', ')))

  x
}


#' Compute VIBER Binomial fits from counts data.
#'
#' @description
#'
#' This function wraps a call to the \code{variational_fit} function
#' from package \code{VIBER}, which implements a variational Binmial
#' mixture model to cluster read counts. See the package manual for
#' \code{VIBER} to see the parameters that can be forwarded, and the
#' plots that can be generated afterwards.
#'
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param samples The samples on which VIBER will be run.
#' @param ... Extra parameters that will be forwarded to
#' \code{variational_fit} function from package \code{VIBER}.
#'
#' @return A new object with a VIBER fit for the required
#' samples in \code{x$fit.Binomial}.
#'
#' @import VIBER
#'
#' @export
#'
#' @examples
#' TODO
mobster_fit_VIBER = function(x, samples = x$samples, ...)
{
  require(VIBER)

  # Prepare inputs
  nv = NV_table(x, samples = samples)
  dp = DP_table(x, samples = samples)

  colnames(nv) = gsub(pattern = '.NV', replacement = '', colnames(nv))
  colnames(dp) = gsub(pattern = '.DP', replacement = '', colnames(dp))

  nv = nv[, colnames(dp)]

  pio::pioTit("Input for variational Binomial clustering with VIBER")

  pio::pioStr("NV values ", '')
  print(nv)
  pio::pioStr("DP values", '')
  print(dp)

  # Wrap fit call
  x$fit.Binomial = VIBER::variational_fit(
    x = nv %>% select(-id),
    y = dp %>% select(-id),
    ...
  )

  # Log update
  x = logOp(x, paste0("Fit VIBER to ", paste0(samples, collapse = ', ')))

  x
}


#' Compute sciClone fits for the available samples.
#'
#' @description
#'
#' This function wraps a call to the sciClone main fitting function
#' from package \code{sciClone}. The function applies the fitting
#' function to the samples stored in the dataset \code{x}, and the
#' returned object contains can be inspected
#' with the functions available in the \code{sciClone} package.
#'
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ... Extra parameters that will be forwarded to
#' \code{mobster_fit} function from package \code{mobster}.
#' @param minimumDepth Mutations below this threshold will be removed.
#' @param maximumClusters Maximum number of clusters to return.
#'
#' @return A new object with sciClone fits in \code{x$fit.sciClone}.
#'
#' @import sciClone
#'
#' @export
#'
#' @examples
#' TODO
mobster_fit_sciClone = function(x,
                                minimumDepth = 0,
                                maximumClusters = 2 * length(x$samples),
                                ...)
{
  inputs = convert_sciClone_input(x)

  names(inputs$CN) = x$samples
  names(inputs$MUTS) = x$samples

  library(sciClone)

  pio::pioHdr("Fitting read counts with the Binomial model available in sciClone",
              toPrint = c(
                `Minimum depth` = minimumDepth,
                `Maximum clusters` = maximumClusters
              ))

  sciClone.fit = sciClone(
    vafs = inputs$MUTS,
    copyNumberCalls = inputs$CN,
    sampleNames = x$samples,
    minimumDepth = minimumDepth,
    maximumClusters = maximumClusters,
    clusterMethod = 'bmm',
    ...
  )

  x$fit.sciClone = sciClone.fit

  # Log update
  x = logOp(x, paste0("Fit sciClone to", x$samples, collapse = ', '))

  x
}

# Converter for sciCLone inputs
convert_sciClone_input = function(x)
{
  seg_ids = unique(x$segments$id)

  # CNA
  CN = NULL
  for (s in x$samples)
  {
    CN.copies = minor(x, seg_id = seg_ids, samples = s) +
      Major(x, seg_id = seg_ids, samples = s)

    CN.segment = x$segments %>%
      select(chr, from, to) %>%
      distinct()
    colnames(CN.segment) = c("chr", 'start', 'stop')

    CN.segment$chr = as.numeric(sapply(CN.segment$chr, function(w)
      substr(w, 4, nchar(w))))

    CN.segment$segment_mean = CN.copies

    CN = append(CN, list(as.data.frame(CN.segment)))
  }

  names(CN) = x$samples

  # MUTAITONS
  locs = x$locations %>% spread(variable, value)

  MUTS = NULL

  for (s in x$samples)
  {
    sample.data = bind_rows(VAF(x, samples = s),
                            DP(x, samples = s),
                            NV(x, samples = s))

    sample.data = sample.data %>% spread(variable, value)
    sample.data = sample.data[complete.cases(sample.data),]

    sample.data = sample.data %>% select(id, DP, NV, VAF)
    sample.data = sample.data %>%
      mutate(refCount = DP - NV,
             varCount = NV,
             vaf = VAF * 100) %>%
      select(id, refCount, varCount, vaf)

    sample.locs = locs %>%
      filter(id %in% sample.data$id) %>%
      mutate(start = from) %>%
      select(id, chr, start)
    sample.locs$start = as.numeric(sample.locs$start)

    df = full_join(sample.data, sample.locs, by = 'id')
    df = df %>% select(chr, start, refCount, varCount, vaf, id)

    df$chr = as.numeric(sapply(df$chr, function(w)
      substr(w, 4, nchar(w))))

    MUTS = append(MUTS,
                  list(as.data.frame(df)))
  }

  list(CN = CN, MUTS = MUTS)
}


# =-=-=-=-=-=-=-=-=-=-=-
# pyClone. Roth et al.
# =-=-=-=-=-=-=-=-=-=-=-


#' Compute pyClone fits for the available samples.
#'
#' @description
#'
#' This function wraps a call to pyClone, a python
#' tool which must be installed and accessyble through a
#' system call within R. The function applies the fitting
#' function to the samples stored in the dataset \code{x}, and the
#' returned object is built scanning pyClone outputs and loading
#' the relevant information.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param iterations Number of pyClone iterations (MCMC).
#' @param burnin Burnin for pyClone MCMC.
#' @param seed Random seed.
#'
#'
#' @return A new object with pyClone fits in \code{x$fit.pyClone}.
#'
#' @export
#'
#' @examples
#' TODO
MOBSTER_fit_pyClone <-
  function(x, iterations=10000, burnin=1000, seed=100) {

    dataset = x

    # Check inputs:
    if (!is.numeric(iterations) | !is.numeric(burnin)){
      stop("Parameters 'interations' or 'burnin' not numeric.\n")
    }

    if (is.na(iterations) | is.na(burnin)) {
      stop("Parameters 'interations' and 'burnin' can not be 'NA'.\n")
    }

    if (!"mbst_data" %in% class(dataset)) {
      stop("Parameters 'dataset' has to be a mobster object.\n")
    }


    # Export data:

    pyclone_workdir <- tempfile(pattern="")
    dir.create(pyclone_workdir, recursive=TRUE, showWarnings=FALSE)

    samples <- dataset$samples

    sfiles <- file.path(pyclone_workdir, paste0(samples, ".tsv"))  # output file names.

    sdata <- lapply(samples, function(sample) { # construct list of data
      mut_ids = keys(dataset)

      nv = NV(dataset, ids=mut_ids, samples=sample)
      nv = nv[match(mut_ids, nv$id),] # ensure correct order

      dp = DP(dataset, ids=mut_ids, samples=sample)
      dp = dp[match(mut_ids, dp$id),] # ensure correct order


      # return
      data.frame(
        mutation_id = mut_ids,
        ref_counts = dp$value - nv$value,
        var_counts = nv$value,
        normal_cn = 2,
        minor_cn = 1,
        major_cn = 1
      )
    })


    for (i in seq(samples)) {# write tables
      write.table(sdata[[i]], sfiles[i], sep="\t", row.names=FALSE, quote=FALSE)
    }


    # Sys call to PyClone:
    cmd = paste0("PyClone run_analysis_pipeline", " \\\n",
                 "   --in_files ", paste(sfiles, collapse = " ")," \\\n",
                 "   --samples ", paste(samples, collapse = " "), " \\\n",
                 "   --num_iters ", iterations, " \\\n",
                 "   --burnin ", burnin, " \\\n",
                 "   --seed ", seed, " \\\n",
                 #"   --max_clusters ", max_clusters, " \\\n",
                 "   --working_dir ", pyclone_workdir)

    cat(crayon::red("Calling PyClone:\n\n"), crayon::blue(cmd))
    log = system(cmd, intern = TRUE)


    # Load PyClone outputs:
    trace_files = list.files(file.path(pyclone_workdir, "trace"), full.names = 1)
    trace_data  = lapply(trace_files, read.delim, stringsAsFactor = 0)
    names(trace_data) = gsub("[.]tsv[.]bz2$", "", basename(trace_files))

    cluster_file = file.path(pyclone_workdir, "tables", "cluster.tsv")
    cluster_data = read.delim(cluster_file, stringsAsFactor = 0)

    loci_file = file.path(pyclone_workdir, "tables", "loci.tsv")
    loci_data = read.delim(loci_file, stringsAsFactor = 0)

    results = list(loci = loci_data, clusters = cluster_data, traces = trace_data)


    # Clean up temp dir:
    unlink(pyclone_workdir, recursive = TRUE)


    # return
    dataset$fit.pyClone = results
    return(dataset)
  }


# # Title
# #
# # @param x
# # @param ...
# #'
# #' @return
# #' @export
# #'
# #' @examples
# mobster_fit_BMM = function(x, ...)
#' {
#'   DP_df = NV_df = NULL
#'
#'   for(s in x$samples)
#'   {
#'     data = DP(x, samples = s) %>%
#'       tidyr::spread(variable, value)
#'
#'     colnames(data)[3] = s
#'     DP_df = bind_cols(DP_df, data[, 3])
#'
#'     data = NV(x, samples = s) %>%
#'       tidyr::spread(variable, value)
#'
#'     colnames(data)[3] = s
#'     NV_df = bind_cols(NV_df, data[, 3])
#'
#'   }
#'
#'   x$fit.BMM = vbdbmm::bmm_mv_em(
#'     successes = as.matrix(NV_df),
#'     trials = as.matrix(DP_df),
#'     ...)
#'
#'   # Log update
#'   x = logOp(x, paste0("Fit BMM to ", x$samples, collapse = ', '))
#'
#'   x
#' }

