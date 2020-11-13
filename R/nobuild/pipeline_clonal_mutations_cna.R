pipeline_clonal_mutations = function(x)
{

}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse tumour sample with MOBSTER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyse_mobster = function(x,
                           cutoff_miscalled_clonal,
                           cna_map,
                           cutoff_lv_assignment,
                           chromosomes,
                           ...)
{
  cli::cli_h1("Analysing tumour sample with MOBSTER")
  cat('\n')

  #  MOBSTER fit of the input data
  mobster_fit_tumour = mobster::mobster_fit(x,
                                            ...)

  output = mobster::Clusters(mobster_fit_tumour$best) %>%
    rename(tumour.cluster = cluster)

  # Determine clonal cluster
  clonal_cluster = guess_mobster_clonal_cluster(
    mobster_fit_tumour = mobster_fit_tumour,
    cutoff_miscalled_clonal = cutoff_miscalled_clonal,
    use_heuristic = analysis_mode(cna_map) == "NO_CNA",
    chromosomes = chromosomes)

  stopifnot(length(clonal_cluster) > 0)

  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
  # Tumour purity: estimated based on the fact that we could map or no the mutations to CNA
  estimated_tumour_purity = NULL
  if(all(is.null(cna_map)))
  {
    # If we could not, then we just assume everything is diploid, therefore
    # therefore 2 * the Beta peak of the clonal cluster

    estimated_tumour_purity = mobster_fit_tumour$best$Clusters %>%
      filter(cluster %in% clonal_cluster, type == 'Mean') %>%
      pull(fit.value) %>% mean * 2

  }
  else
  {
    # Otherwise we do something a little bit smarter, which is normalising for segment's ploidy.

    # We take the mean of the clonal cluster
    ccluster_mean = mobster_fit_tumour$best$Clusters %>%
      filter(cluster %in% clonal_cluster, type == 'Mean') %>%
      pull(fit.value) %>%
      mean

    # We found the karyotype of the mutations that we used (we take the first one, because they are all the same)
    used_karyotype = strsplit(x$karyotype[1], ':')[[1]]

    # By design, we should have correctly taken as coonal cluster the set of mutations that happened before
    # aneuploidy. Also, again by the fact that we use only simple karyotypes, we assume that the mutation is
    # present in a number of copies that match the actual Major allele (M). This is the same thing as saying
    # that the mutations happened *before* aneuploidy
    estimated_tumour_purity = CNAqc:::purity_estimation_fun(
      v = ccluster_mean, # Observed VAF
      m = used_karyotype[2] %>% as.numeric,
      M = used_karyotype[1] %>% as.numeric,
      mut.allele = 2
    )
  }
  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

  # List of clonal mutations in the tumour, with LV > cutoff_lv_assignment
  clonal_tumour = NULL
  repeat
  {
    clonal_tumour = mobster::Clusters(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) %>%
      filter(cluster %in% clonal_cluster) %>%
      pull(id)

    if(length(clonal_tumour) > 20 | cutoff_lv_assignment < 0) break

    cutoff_lv_assignment = cutoff_lv_assignment - 0.03

    # cat("Dynamic adjustment: ", cutoff_lv_assignment, " n = ", length(clonal_tumour))
  }

  # Plot
  figure = plot_mobster_fit(mobster_fit_tumour, cutoff_lv_assignment, clonal_cluster)

  return(
    list(
      output = output,
      fit = mobster_fit_tumour,
      plot = figure,
      clonal_cluster = clonal_cluster,
      clonal_mutations = clonal_tumour,
      estimated_purity = estimated_tumour_purity,
      cutoff_lv_assignment = cutoff_lv_assignment
    )
  )
}
