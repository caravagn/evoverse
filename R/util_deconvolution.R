deconvolution_mobster_karyotypes = function(mutations,
                                            karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                            min_muts = 50,
                                            ...)
{

  timeable = karyotypes

  # Timeable ones, processed one by one with MOBSTER
  # VAF required in 0/1
  mfits = lapply(
    timeable,
    function(k)
    {
      cat("\n")
      cli::cli_process_start("MOBSTER clustering of mutations with karyotype {.field {k}}")
      cat("\n")

      kmuts = mutations %>%
        filter(karyotype == k, VAF > 0, VAF < 1)

      if(nrow(kmuts) < min_muts) {
        cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
        cli::cli_process_failed()

        return(NULL)
      }

      mfit = mobster::mobster_fit(kmuts, ...)

      cat("\n")
      cli::cli_process_done()

      return(mfit)
    }
  )
  names(mfits) = timeable

  mfits
}

deconvolution_mobster_CCF = function(cna_obj,
                                     CCF_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                     min_muts = 50,
                                     ...)
{

  if(all(is.null(CCF_karyotypes))) stop("CCF_karyotypes should not be null")

  cna_obj = CNAqc::compute_CCF(cna_obj, karyotypes = CCF_karyotypes)

  # Timeable ones, processed one by one with MOBSTER
  # VAF required in 0/1
  mfits = lapply(
    timeable,
    function(k)
    {
      cat("\n")
      cli::cli_process_start("MOBSTER clustering of mutations with karyotype {.field {k}}")
      cat("\n")

      kmuts = mutations %>%
        filter(karyotype == k, VAF > 0, VAF < 1)

      if(nrow(kmuts) < min_muts) {
        cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
        cli::cli_process_failed()

        return(NULL)
      }

      mfit = mobster::mobster_fit(kmuts, ...)

      cat("\n")
      cli::cli_process_done()

      return(mfit)
    }
  )
  names(mfits) = timeable

  mfits
}

qc_deconvolution_mobster = function(x, type)
{
  qc_model = NULL
  input = mobster::to_string(x)

  # Specialise QC and data
  if(type == 'T') {
    qc_model = evoverse::qc_timing_model

    x$QC_type = 'timing_logistic_GEL_CRC'

    # Fail monoclonal for timing!
    if(x$Kbeta == 1)
    {
      x$QC <- 'FAIL'
      x$QC_prob <- 0

      return(x)
    }

    input = input  %>%
      mutate(ratioM = abs((Mean_C1 / Mean_C2) - 2),
             ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
             ratioN = (max(N_C1, N_C2) +1) / (min(N_C2, N_C1) + 1),
             minpi = min(pi_C1, pi_C2))
  }

  # Specialise QC and data
  if(type == 'D')
  {
      if(x$Kbeta == 1) qc_model = evoverse::qc_deconvolution_monoclonal
      else qc_model = evoverse::qc_deconvolution_polyclonal

       x$QC_type = 'timing_logistic_GEL_CRC'

      input = input %>%
        select(starts_with('sse'),
               Variance_C1,
               Variance_C2,
               reduced.entropy,
               entropy)
  }

  if(is.null(qc_model))  stop("QC available: (T) Timing, (D) Deconvolution.")

  # Make predictions
  x$QC_prob <- qc_model %>% predict(input, type = "response")
  x$QC <- ifelse(x$QC_prob > 0.5, "PASS", "FAIL")

  return(x)
}

qc_mobster_plot = function(x)
{
  if(all(is.null(x))) return(ggplot() + geom_blank())

  qc = ifelse(x$QC == "FAIL", "indianred3", 'forestgreen')

  mobster::plot.dbpmm(x) +
    labs(title = bquote("AUTO QC " ~ .(x$QC) ~ p["PASS"] ~
                          '=' ~ .(x$QC_prob))) +
    theme(title = element_text(color = qc),
          panel.border = element_rect(
            colour = qc,
            fill = NA,
            size = 5
          ))
}
