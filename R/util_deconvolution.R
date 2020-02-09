######################################################
### INPUT FUNCTIONS
######################################################

# Transforms input data into a CNAqc object, if required it downsamples the data.
deconvolution_prepare_input = function(mutations, cna, purity, N_max, min_VAF)
{
  # x must be a dataset for MOBSTER
  if(!is.data.frame(mutations) & !is.matrix(mutations))
    stop("mutations: Must use a dataset for MOBSTER!")

  if(!is.null(cna) & !is.data.frame(cna) & !is.matrix(cna))
    stop("cna: Must use a dataset for MOBSTER!")

  if(!is.null(purity) & (purity > 1 | purity <= 0))
    stop("Purity must be in 0/1.")

  # VAF above cutoff
  if(any(mutations$VAF < min_VAF, na.rm = T))
  {
    cli::cli_alert_info("Removing {.value {sum(mutations$VAF < min_VAF, na.rm = T)}} mutations with VAF < {.value {min_VAF}}.")
    mutations = mutations %>% dplyr::filter(VAF > min_VAF)
  }

  # Downsample data
  if(mutations %>% nrow > N_max) {

    cli::boxx(paste0("n = ", mutations %>% nrow,  " mutations. n > ", N_max, " , downsampling input.")) %>% cat
    cat('\n')

    mutations = mutations %>% dplyr::sample_n(N_max)
  }

  # Apply CNA mapping and retain only mappable mutations
  return(CNAqc::init(mutations, cna, purity))
}

######################################################
### DECONVOLUTION WITH MOBSTER WOTH AUTO QC AND AUTO-RERUN
######################################################

# # Fit every subset of mutations (defined by karyotypes) with MOBSTER. Skip those with < min_muts.
smartrun_mobster_qc = function(x, QC_type, model_description = "My model", ...)
{
  # First, attempt a mobster fit with the required parameters
  mfit = mobster::mobster_fit(x, description = model_description, ...)$best

  # Then, QC the model fit
  mfit = evoverse:::qc_deconvolution_mobster(mfit, type = QC_type)

  cat("\n")
  if(mfit$QC == "PASS") cli::cli_alert_success("{.field {green(model_description)}}: QC PASS. p = {.value {mfit$QC_prob}}")
  else cli::cli_alert_danger("{.field {red(model_description)}}: QC FAIL. p = {.value {mfit$QC_prob}}")
  cat("\n")

  if(is.null(QC_type) |  mfit$QC == "PASS") return(mfit)

  cat("\n")
  cat(
    cli::boxx(
      paste0("Auto-QC (type '", QC_type, "') failed, attempting another single run with default parameters"),
      padding = 1,
      float = 'center',
      background_col = "brown")
  )
  cat("\n")

  # Second attempt
  mfit = mobster::mobster_fit(x, description = model_description, parallel = FALSE)

  # Then, QC the model fit
  mfit = evoverse:::qc_deconvolution_mobster(mfit$best, type = QC_type)

  cat("\n")
  if(mfit$QC == "PASS") cli::cli_alert_success("{.field {green(model_description)}}: QC PASS. p = {.value {mfit$QC_prob}}")
  else cli::cli_alert_danger("{.field {red(model_description)}}: QC FAIL. p = {.value {mfit$QC_prob}}")
  cat("\n")

  return(mfit)
}

######################################################
### DECONVOLUTION WITH MOBSTER/ BMIX OF VAF VALUES
######################################################


# Fit every subset of mutations (defined by karyotypes) with MOBSTER. Skip those with < min_muts.
deconvolution_mobster_karyotypes_VAF = function(mutations,
                                            karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                            min_muts = 50,
                                            QC_type = NULL,
                                            ...)
{
  # Timeable ones, processed one by one with MOBSTER
  # VAF required in 0/1
  mfits = lapply(karyotypes,
                 function(k)
                 {
                   cat("\n")
                   cli::cli_h1("MOBSTER clustering mutations with karyotype {.field {k}}")
                   cat("\n")

                   kmuts = mutations %>% filter(karyotype == k, VAF > 0, VAF < 1)

                   if (nrow(kmuts) < min_muts) {
                     cli::cli_alert_warning("Less than {.value {min_muts}} mutations, will no fit.")

                     return(NULL)
                   }

                   mfit = smartrun_mobster_qc(kmuts, QC_type, model_description = paste0("Raw VAF in ", k),  ...)

                   return(mfit)
                 })
  names(mfits) = karyotypes

  mfits
}

deconvolution_nonread_counts_bmix_VAF = function(x,
                                                 karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                                 min_muts = 50)
{
  runner = function(x, k)
  {
    if(all(is.null(x[[k]]))) {
      cli::cli_alert_warning("Requested to fit read counts from {.field {k}}, but the required MOBSTER fit is null.")
      # Some error message?
      return(NULL)
    }

    cat("\n")
    cli::cli_h1("BMix clustering non-tail mutations with karyotype {.field {k}}")
    cat("\n")

    # Get non_tail mutations
    non_tail = mobster::Clusters(x[[k]]) %>%
      dplyr::filter(cluster != "Tail") %>%
      data.frame

    # Not enough non-tail mutations
    if(nrow(non_tail) < min_muts) {
      cli::cli_alert_warning("Requested to fit read counts from {.field {k}}, but therere are less than {.value {min_muts}} available non-tail mutations.")

      return(NULL)
    }

    # Fit: use 2 times the clusters detected at max

    Kbeta = min(x[[k]]$Kbeta * 2, 4)

    fit_readcounts = BMix::bmixfit(
      non_tail %>% dplyr::select(NV, DP),
      K.BetaBinomials = 0,
      K.Binomials = 1:Kbeta,
      samples = 3)

    fit_readcounts$input = non_tail
    fit_readcounts
  }

  lapply(karyotypes, runner, x = x)
}


######################################################
### DECONVOLUTION WITH MOBSTER/ BMIX OF CCF VALUES
######################################################

deconvolution_mobster_CCF = function(cna_obj,
                                     CCF_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                     min_muts = 50,
                                     min_VAF = 0.05,
                                     QC_type = NULL,
                                     ...)
{
  if (all(is.null(CCF_karyotypes)) | length(CCF_karyotypes) == 0) {
    warning("CCF_karyotypes are null or empty")
    return(NULL)
  }

  cat("\n")
  cli::cli_h1("MOBSTER clustering CCF for mutations with karyotype(s) {.field {CCF_karyotypes}}")
  cat("\n")

  # print(summary(mutations))
  available_karyo = cna_obj$n_karyotype[CCF_karyotypes]
  available_karyo = available_karyo[!is.na(available_karyo)]

  if(!all(CCF_karyotypes %in% names(available_karyo)))
  {
    CCF_karyotypes = intersect(CCF_karyotypes, names(available_karyo))

    cli::cli_alert_info("Some karyotypes are not avilable, using {.field {CCF_karyotypes}}")
  }

  cna_obj = CNAqc::compute_CCF(cna_obj, karyotypes = CCF_karyotypes)
  mutations = Reduce(bind_rows,
                     lapply(cna_obj$CCF_estimates, function(x)
                       x$mutations))

  # We remove things that make no sense.. NA for CCF
  mutations = mutations %>%
    dplyr::filter(!is.na(CCF))

  if (nrow(mutations) < min_muts)
  {
    cli::cli_alert_warning("Less than {.value {min_muts}} mutations, will not fit CCF data.")
    return(list(fits = NULL, cna_obj = cna_obj))
  }

  # Scale CCF - first remove anything which exceeds 2
  if(any((mutations$CCF/2) > 1, na.rm = T))
  {
    cli::boxx(
      paste0(
        "n = ", sum(mutations$CCF/2 > 1),
        " mutation(s) with CCF/2 > 1 will be removed (50% adjusted VAF).")) %>%
      cat
    cat('\n')
  }

  # Then we do something else to reflect the fact that the VAF
  # and CCF values are related, and if any filter is applied to VAF data,
  # then the same filter should be applied to CCF
  M_VAF = mutations %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(m = min(VAF, na.rm = T)) %>%
    tidyr::separate(karyotype, into = c("A", "B"), sep =':') %>%
    dplyr::mutate(
      min_CCF = CNAqc:::ccf_adjustment_fun(min_VAF, as.numeric(A), as.numeric(B), cna_obj$purity, mut.allele = 1),
      min_VAF = min_VAF
    )

  M = M_VAF %>%
    dplyr::summarise(M = max(min_CCF, na.rm = T)) %>%
    dplyr::pull()
  M = M/2

  cat('\n')
  cli::boxx(paste0("Minimum adjusted VAF (=CCF/2) value ", round(M, 4), ", from VAF/CCF profiling")) %>% cat
  cat('\n')

  # CCF values divided by 2 to get adjusted VAF
  mutations = mutations %>%
    dplyr::mutate(
      VAF_raw = VAF,
      VAF = CCF / 2) %>%
    dplyr::filter(VAF > M, VAF < 1)

  if (nrow(mutations) < min_muts)
  {
    cli::cli_alert_warning("Less than {.value {min_muts}} mutations, will not fit CCF data.")
    return(list(fits = NULL, cna_obj = cna_obj))
  }

  return(list(
    fits = smartrun_mobster_qc(mutations, QC_type, model_description = "CCF values", ...),
    cna_obj = cna_obj
  ))
}

deconvolution_nonread_counts_bmix_CCF = function(x, min_muts)
{
  if(all(is.null(x[["CCF"]]))) return(NULL)

  x = x[["CCF"]]

  cat("\n")
  cli::cli_h1("BMix clustering non-tail mutations from CCF data")
  cat("\n")

  # Get non_tail mutations
  non_tail = mobster::Clusters(x) %>%
      dplyr::filter(cluster != "Tail") %>%
      data.frame

  # Not enough non-tail mutations
  if(nrow(non_tail) < min_muts) {
    cli::cli_alert_warning("Requested to fit read counts from {.field {k}}, but therere are less than {.value {min_muts}} available non-tail mutations.")

    return(list(`CCF` = NULL))
  }

  # CCF have been computed before, and adjusted VAF values (CCF/2) have been
  # used to cluster with MOBSTER. Now we hold fixed the sequencing depth DP,
  # and derive from NV/DP=VAF that NV=VAF*DP, and adjust NV and NR constrained to NV+NR=DP
  # (no coverage change, trials).
  non_tail = mobster::Clusters(x) %>%
    dplyr::filter(cluster != "Tail") %>%
    dplyr::mutate(
      NV = floor(VAF * DP),
      NR = DP - NV,
      NV = ifelse(NV < 0, 0, NV),
      NV = ifelse(NV > DP, DP, NV),
      NR = ifelse(NR < 0, 0, NR),
      NR = ifelse(NR > DP, DP, NR),
      VAF = NV / (NV + NR)
    ) %>%
    data.frame

  Kbeta = min(x$Kbeta * 2, 4)

  fit_readcounts = BMix::bmixfit(
    non_tail %>% dplyr::select(NV, DP),
    K.BetaBinomials = 0,
    K.Binomials = 1:Kbeta,
    samples = 3)

  fit_readcounts$input = non_tail
    fit_readcounts
}


# Perform QC of a MOBSTER fit, either based on the Timing or Deconvolution classfiers
qc_deconvolution_mobster = function(x, type)
{
  qc_model = NULL
  if(all(is.null(x))) return(NULL)

  input = mobster::to_string(x)

  # Specialise QC and data
  if (type == 'T') {
    qc_model = evoverse::qc_timing_model

    x$QC_type = 'timing_logistic_GEL_CRC'

    # Fail monoclonal for timing!
    if (x$Kbeta == 1)
    {
      x$QC <- 'FAIL'
      x$QC_prob <- 0

      return(x)
    }

    input = input  %>%
      dplyr::mutate(
        ratioM = abs((Mean_C1 / Mean_C2) - 2),
        ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
        ratioN = (max(N_C1, N_C2) + 1) / (min(N_C2, N_C1) + 1),
        minpi = min(pi_C1, pi_C2)
      )
  }

  # Specialise QC and data
  if (type == 'D')
  {
    if (x$Kbeta == 1) {
      qc_model = evoverse::qc_deconvolution_monoclonal

      input = input %>%
        dplyr::select(starts_with('sse'),
               Variance_C1,
               reduced.entropy,
               entropy,
               tail)

    }
    else {
      qc_model = evoverse::qc_deconvolution_polyclonal

      input = input %>%
        dplyr::select(starts_with('sse'),
               Variance_C1,
               Variance_C2,
               reduced.entropy,
               tail,
               entropy)
    }

    x$QC_type = 'deconvolution_logistic_Hartwig'
  }

  if (is.null(qc_model))
    stop("QC available: (T) Timing, (D) Deconvolution.")

  # Make predictions
  x$QC_prob <- qc_model %>% predict(input, type = "response") %>% round(digits = 3)
  x$QC <- ifelse(x$QC_prob > 0.5, "PASS", "FAIL")

  return(x)
}



qc_mobster_plot = function(x)
{
  if (all(is.null(x)))
    return(ggplot() + geom_blank())

  qc = ifelse(x$QC == "FAIL", "indianred3", 'forestgreen')

  mobster::plot.dbpmm(x) +
    labs(title = bquote("QC " ~ .(x$QC) ~ p["PASS"] ~
                          '=' ~ .(x$QC_prob))) +
    theme(title = element_text(color = qc),
          panel.border = element_rect(
            colour = qc,
            fill = NA,
            size = 5
          ))
}


# Assemble the plot
deconvolution_plot_assembly = function(mobster_fits, cna_obj, bmix_fits, figure_caption, figure_title)
{
  groups = names(mobster_fits)
  empty_panel = ggplot() + geom_blank()

  # MOBSTER plots, sourrounded by a coloured box by QC
  mob_fits_plot = lapply(groups, function(y) evoverse:::qc_mobster_plot(mobster_fits[[y]]))

  # CNA plot
  cna_plot = empty_panel
  if(!is.null(cna_obj)) cna_plot = CNAqc::plot_segments(cna_obj, circular = TRUE)

  # Top panel: CNA + MOBSTER
  mob_fits_plot =  append(list(cna_plot), mob_fits_plot)
  figure = ggpubr::ggarrange(
    plotlist = mob_fits_plot,
    nrow = 1,
    ncol = length(groups) + 1,
    labels = c("CNA", groups)
    )

  # If there is a second panel, we put it below
  if(!all(is.null(bmix_fits)))
  {
    # BMIx: a panel like the one above, same dimension
    bmix_panel = lapply(
      groups,
      function(x)
      {
        if (all(is.null(bmix_fits[[x]]))) return(empty_panel)
        BMix::plot_clusters(bmix_fits[[x]], bmix_fits[[x]]$input %>% dplyr::select(NV, DP))
      })

    bmix_panel = ggarrange(plotlist = append(list(empty_panel), bmix_panel),
                           nrow = 1,
                           ncol = length(groups) + 1,
                           labels = c("", groups))

    # Assemble a one-page final figure with both MOBSTER and BMix panels
    figure =  ggpubr::ggarrange(
      figure,
      bmix_panel,
      nrow = 2,
      ncol = 1
    )
  }

  # Set the figure title and captions
  figure = ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(bquote(bold("Dataset. ") ~ .(figure_title)), hjust = 0, x = 0, size = 15),
    bottom = ggpubr::text_grob(bquote(.(figure_caption)), hjust = 0, x = 0, size = 8)
  )

  return(figure)
}


# Extract a table of mutations assignments: MOBSTER + BMix fits
deconvolution_table_assignments = function(mobster_fits, bmix_fits)
{
  groups = names(mobster_fits)

  # First, gather all VAF-based clusters
  summary_table = Reduce(bind_rows,
                         lapply(groups,
                                function(x)
                                {
                                  mfit = mobster_fits[[x]]
                                  bfit = bmix_fits[[x]]

                                  # If there is no fit -> empty table
                                  if (all(is.null(mfit)))
                                    return(NULL)

                                  # Otherwise take the VAF clustering for this karyotype
                                  mobster_tab = mobster::Clusters(mfit) %>% dplyr::mutate(karyotype = x)

                                  if(all(is.null(bfit))) return(mobster_tab)

                                  bmix_tab = BMix::Clusters(bfit, bfit$input) %>%
                                    dplyr::select(chr, from, to, ref, alt, karyotype, cluster)

                                  mobster_tab %>%
                                    dplyr::full_join(bmix_tab,
                                                     by = c('chr', 'from', 'to', 'ref', 'alt', 'karyotype')) %>%
                                    dplyr::rename(cluster = cluster.x, BMix_cluster = cluster.y) %>%
                                    dplyr::select(-ends_with('.y'))
                                }))

  # # If there is CCF data out, take genome coordinates, merge the CCF value and cluster
  # if ("CCF" %in% groups)
  # {
  #   mCCF_output = mobster::Clusters(mobster_fits$CCF) %>%
  #     dplyr::select(chr, from, to, ref, alt, CCF, cluster) %>%
  #     dplyr::rename(CCF_cluster = cluster)
  #
  #   summary_table = summary_table %>%
  #     dplyr::full_join(mCCF_output, by = c('chr', 'from', 'to', 'ref', 'alt'))
  #
  #   bCCF_output = BMix::Clusters(bmix_fits$CCF, bmix_fits$CCF$input) %>%
  #     dplyr::select(chr, from, to, ref, alt, cluster) %>%
  #     dplyr::rename(CCF_BMix_cluster = cluster)
  #
  #   summary_table = summary_table %>%
  #     dplyr::full_join(bCCF_output, by = c('chr', 'from', 'to', 'ref', 'alt'))
  # }

  return(summary_table)
}

# Extract a table of summary stats: MOBSTER + BMix fits
deconvolution_table_summary = function(mobster_fits, bmix_fits)
{
  groups = names(mobster_fits)

  summary_table = lapply(
    groups,
    function(k)
    {
      # No Mobster, no table
      if(mobster_fits[[k]] %>% is.null %>% all) return(NULL)

      x = mobster_fits[[k]]

      # Mobster summary, plus QC results
      m_tab = bind_cols(
        mobster::to_string(x),
        data.frame(QC = x$QC, QC_prob = x$QC_prob, QC_type = x$QC_type, karyotype = k, stringsAsFactors = F)
      )

      if(all(is.null(bmix_fits[[k]]))) return(m_tab)

      # BMix summary
      b_tab = BMix::to_string(bmix_fits[[k]])
      colnames(b_tab) = paste0('BMix_', colnames(b_tab))

      return(cbind(m_tab, b_tab) %>% as_tibble())
      })

  return(Reduce(bind_rows, summary_table))
}


# bckg_color = RColorBrewer::brewer.pal(
#   n = length(atab$CCF_cluster %>% unique), name = 'Set1') %>%
#   alpha(.8)
#
# atab %>%
#   ggplot(
#     aes(x = CCF, fill = CCF_cluster)
#   ) +
#   geom_histogram(binwidth = 0.01, size = .1) +
#   scale_fill_manual(values = bckg_color) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_col
#   mobster:::my_ggplot_theme() +
#   facet_wrap(~cluster, ncol = 1)
#   #acet_grid(cluster ~ CCF_cluster)


