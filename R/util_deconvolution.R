deconvolution_prepare_input = function(mutations, cna, purity, N_max)
{
  # x must be a dataset for MOBSTER
  if(!is.data.frame(mutations) & !is.matrix(mutations))
    stop("mutations: Must use a dataset for MOBSTER!")

  if(!is.null(cna) & !is.data.frame(cna) & !is.matrix(cna))
    stop("cna: Must use a dataset for MOBSTER!")

  if(!is.null(purity) & (purity > 1 | purity <= 0))
    stop("Purity must be in 0/1.")

  # Downsample data
  if(mutations %>% nrow > N_max) {

    cli::boxx(paste0("n = ", mutations %>% nrow,  " mutations. n > ", N_max, " , downsampling input.")) %>% cat
    cat('\n')

    mutations = mutations %>% dplyr::sample_n(N_max)
  }

  # Apply CNA mapping and retain only mappable mutations
  cna_obj = NULL
  if(!is.null(cna))
  {
    cli::cli_h1("Found CNA calls, retaining mutations mapping to available segments")
    cat("\n")

    cna_obj = CNAqc::init(mutations, cna, purity)
    mutations = cna_obj$snvs
  }

  return(list(mutations = mutations, cna = cna, purity = purity, cna_obj = cna_obj))
}

deconvolution_mobster_karyotypes = function(mutations,
                                            karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                            min_muts = 50,
                                            ...)
{
  timeable = karyotypes


  # Special case -- if karyotypes is null all data are used
  if (all(is.null(timeable)))
  {
    k = 'wg'

    cat("\n")
    cli::cli_h1("MOBSTER clustering all mutations (no subsetting by karyotype)")
    cat("\n")

    kmuts = mutations %>% filter(VAF > 0, VAF < 1)

    if (nrow(kmuts) < min_muts) {
      cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
      # cli::cli_process_failed()

      return(NULL)
    }

    mfit = mobster::mobster_fit(kmuts, ...)

    return(list(`wg` = mfit))
  }


  # Timeable ones, processed one by one with MOBSTER
  # VAF required in 0/1
  mfits = lapply(timeable,
                 function(k)
                 {
                   cat("\n")
                   cli::cli_h1("MOBSTER clustering mutations with karyotype {.field {k}}")
                   cat("\n")

                   kmuts = mutations %>%
                     filter(karyotype == k, VAF > 0, VAF < 1)

                   if (nrow(kmuts) < min_muts) {
                     cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
                     # cli::cli_process_failed()

                     return(NULL)
                   }

                   mfit = mobster::mobster_fit(kmuts, ...)

                   # cat("\n")
                   # cli::cli_process_done()

                   return(mfit)
                 })
  names(mfits) = timeable

  mfits
}

deconvolution_mobster_CCF = function(cna_obj,
                                     CCF_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                     min_muts = 50,
                                     ...)
{
  if (all(is.null(CCF_karyotypes)))
    stop("CCF_karyotypes should not be null")

  # cli::cli_process_start("Performing deconvolution with CCF from karyotypes {.field {CCF_karyotypes}}")

  cat("\n")
  cli::cli_h1("MOBSTER clustering CCF for mutations with karyotype(s) {.field {CCF_karyotypes}}")
  cat("\n")

  if (!("CCF" %in% colnames(cna_obj$snvs)))
  {
    available_karyo = cna_obj$n_karyotype[CCF_karyotypes]
    available_karyo = available_karyo[!is.na(available_karyo)]

    cna_obj = CNAqc::compute_CCF(cna_obj, karyotypes = names(available_karyo))
    mutations = Reduce(bind_rows,
                       lapply(cna_obj$CCF_estimates, function(x)
                         x$mutations))
  }
  else
    cli::cli_alert_warning("Using CCF annotation already available in the data")

  if(any(mutations$CCF/2 > 1)) {
    cli::boxx(paste0("n = ", sum(mutations$CCF/2 > 1),  " mutations with CCF/2 > will be removed.")) %>% cat
    cat('\n')

  }

  mutations = mutations %>%
    mutate(VAF_raw = VAF,
           VAF = CCF / 2) %>%
    filter(VAF > 0, VAF < 1)

  if (nrow(mutations) < min_muts) {
    cli::cli_alert_warning("Less than {.value {min_muts}} mutations, skipping.")
    # cli::cli_process_failed()

    return(list(fits = NULL, cna_obj = cna_obj))
  }

  return(list(
    fits = mobster::mobster_fit(mutations, ...),
    cna_obj = cna_obj
  ))
}

# Perform QC of a MOBSTER fit
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
                         lapply(groups[groups != 'CCF'],
                                function(x)
                                {
                                  # If there is no fit -> empty table
                                  if (all(is.null(mobster_fits[[x]])))
                                    return(NULL)

                                  # Otherwise take the VAF clustering for this karyotype
                                  mobster_tab = mobster::Clusters(mobster_fits[[x]]) %>%
                                    dplyr::mutate(karyotype = x)

                                  if(all(is.null(bmix_fits[[x]]))) return(mobster_tab)

                                  bmix_tab = BMix::Clusters(bmix_fits[[x]], bmix_fits[[x]]$input) %>%
                                    dplyr::select(chr, from, to, ref, alt, karyotype, cluster)

                                  mobster_tab %>%
                                    dplyr::full_join(bmix_tab,
                                                     by = c('chr', 'from', 'to', 'ref', 'alt', 'karyotype')) %>%
                                    dplyr::rename(cluster = cluster.x, BMix_cluster = cluster.y) %>%
                                    dplyr::select(-ends_with('.y'))
                                }))

  # If there is CCF data out, take genome coordinates, merge the CCF value and cluster
  if ("CCF" %in% groups)
  {
    mCCF_output = mobster::Clusters(mobster_fits$CCF) %>%
      dplyr::select(chr, from, to, ref, alt, CCF, cluster) %>%
      dplyr::rename(CCF_cluster = cluster)

    summary_table = summary_table %>%
      dplyr::full_join(mCCF_output, by = c('chr', 'from', 'to', 'ref', 'alt'))

    bCCF_output = BMix::Clusters(bmix_fits$CCF, bmix_fits$CCF$input) %>%
      dplyr::select(chr, from, to, ref, alt, cluster) %>%
      dplyr::rename(CCF_BMix_cluster = cluster)

    summary_table = summary_table %>%
      dplyr::full_join(bCCF_output, by = c('chr', 'from', 'to', 'ref', 'alt'))
  }

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


