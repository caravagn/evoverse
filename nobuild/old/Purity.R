#' Estimate the purity of
#'
#' @description
#'
#' From the peaks of the VAF distributions
#' for each sample in \code{x}, this function
#' estimates the purity looking for the
#' highest peak in a range \code{peak.range}
#' which is by default (0.2, 0.5, i.e., from 40%
#' to 100% tumour content). An extra parameter
#' \code{m.peak} defines the contour of each
#' peak; change this parameter to get different
#' estimates.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param peak.range A range of VAF values where
#' the peaks are searched for.
#' @param m.peak The neighbour which are compared
#' against when searching for a peak.
#'
#' @return A \code{tibble} with the peaks per sample
#'
#' @export
#'
#' @examples
#' TODO
estimate_purity = function(x, peak.range = c(0.2, 0.5), m.peak = 10)
{
  pio::pioTit(paste("Guessing purity via peak detection in", peak.range[1], '--', peak.range[2], 'for diploid SNVs'))

  message("Warning. This thing makes more sense if your using mostly diploid SNVs \n\n")

  results = NULL

  for(s in x$samples)
  {
    v = VAF(x, samples = s) %>% pull(value)

    # Detect peaks
    v = v[v > peak.range[1]]
    v = v[v < peak.range[2]]

    h = hist(v, breaks = seq(0, 1, 0.01), plot = FALSE)

    peaks = .find_peaks(h$density, m.peak)
    x.peaks = (peaks * 0.01)

    # diploid tumour and normal
    guess = data.frame(samples = s, purity = 2 * x.peaks)

    results = rbind(results, guess)
  }

  results = as_tibble(results)

  results
}


# Peaks detection function
.find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
