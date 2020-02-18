#' Title
#'
#' @param x
#' @param cutoff_lv
#' @param N_min
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' Take one dataset with genomic locations
#' data('PD4120a_breast_sample', package = 'mobster')
#'
#' # Tune the pipeline (fast)
x = pipeline_clonal_mutations(
  PD4120a_breast_sample$best$data,
  cutoff_lv = .7,
  auto_setup = 'FAST', # These go to MOBSTER mobster_fit
  K = 3)

print(x)

pipeline_clonal_mutations = function(x, cutoff_lv = .7, N_min = 20, ...)
{
  # Required input
  required_columns = c('chr', 'from', 'to')
  stopifnot(required_columns %in% colnames(x))

  # This is the default analysis, it
  used_chromsomes = CNAqc::chr_coordinates_GRCh38$chr %>% unique

  #  MOBSTER fit of the input data
  cli::cli_h1("Determining clonal mutations from VAF data (all karyotypes)")
  cat('\n')
  cli::cli_alert_info("Working without CNA data: MOBSTER will use a whole-genome reference")
  cat('\n')

  mobster_fit_tumour = mobster::mobster_fit(x, ...)

  best = mobster_fit_tumour$best

  # Create a table with what seems reasonable
  # 1) high-confidence assignments
  LV_clusters_table = heuristic_mobster_cluster_assignments_dynamic_lv(best, cutoff = cutoff_lv, N_min = N_min)

  # Render the parameter the vector based on dynamic scoring
  cutoff_lv = pio:::nmfy(LV_clusters_table$cluster, LV_clusters_table$cutoff)

  # 2) enrichment likelihood test ~ Are >60% of mutations clustered in <20% of the considered genome (full genome here)
  LOCATION_clusters_table = heuristic_mobster_clustering_enriched_for_chromsome_locations(
    best,
    cutoff_lv,
    chromosomes = used_chromsomes)

  cat('\n')
  cli::cli_alert_warning("Mutations spread across chromsomes (whole-genome)")

  print(
    LOCATION_clusters_table %>% dplyr::select(cluster, prob, p_20p, CLUSTERED)
  )

  # 3) estimate tumour purity. Without mapping mutations to CNA, we assume everything is diploid
  cli::cli_alert_warning("Without mapping mutations to CNA, a 'diploid' adjustment is assumed: purity = 2 * clonal_VAF")
  PURITY_clusters_table = heuristic_mobster_purity_estimation_nocna(best)

  # Select the best clonal cluster
  fit = list()
  fit$fit = best

  # All table data as ouptut joint
  fit$table = LV_clusters_table %>%
    dplyr::full_join(LOCATION_clusters_table, by = 'cluster') %>%
    dplyr::full_join(PURITY_clusters_table, by = 'cluster')

  # Cluster means (to scan in order and determie priority by location)
  cluster_Binomial_means = x$Clusters %>% dplyr::filter(type == 'Mean') %>% dplyr::select(-init.value)

  candidates = fit$table %>%
    dplyr::full_join(cluster_Binomial_means, by = 'cluster') %>%
    dplyr::arrange(desc(fit.value)) %>%
    dplyr::filter(!CLUSTERED)

  # Special case, all seem wrong, jus say and assume 1 is correct
  if(nrow(candidates) == 0)
  {
    candidates = fit$table %>%
      dplyr::full_join(cluster_Binomial_means, by = 'cluster') %>%
      dplyr::arrange(desc(fit.value))

    cli::cli_alert_danger("The mutations in all MOBSTER clusters map to few chromsomes, try to use also CNA data.")
  }

  fit$clonal_cluster = candidates$cluster[1]
  fit$purity =  candidates$purity[1]
  fit$n_clonal = table(x$data$cluster)[fit$clonal_cluster]

  cli::cli_alert_success("Clonal cluster: {.field {fit$clonal_cluster}} with n = {.value {fit$n_clonal}} mutations.")

  if(fit$purity > 1)
    cli::cli_alert_danger("Tumor purity: {.value {fit$purity}} (>1), will be coerced to 1.")

  if(fit$clonal_cluster == "Tail")
    cli::cli_alert_danger("Clonal cluster as 'Tail', does not make sense.")

  return(fit)
}



# x ~ MOBSTER fit
# cutoff ~ desired starting value to determine high-confidence hard clustering assignments
# N_min ~ minimum number of mutations to be considered sufficient for a cluster
heuristic_mobster_cluster_assignments_dynamic_lv = function(x, cutoff, N_min = 20)
{
  stopifnot(cutoff >= 0 & cutoff <= 1)
  stopifnot(N_min >= 0)
  stopifnot(inherits(x, 'dbpmm'))

  # Solves the problem for a single cluster
  dynamic = function(cluster)
  {
    mutations = NULL
    # repeat untill we have the best we can get
    repeat
    {
      # Get top as of `cutoff`
      mutations = mobster::Clusters(x, cutoff_assignment = cutoff) %>%
        dplyr::filter(cluster %in% !!cluster)

      # Stop if we have enough mutations or we have done the best we can
      if(nrow(mutations) > N_min | cutoff < 0) break

      cutoff = cutoff - 0.03
    }

    return(data.frame(cluster = cluster, cutoff = cutoff, N = nrow(mutations), stringsAsFactors = FALSE))
  }

  tibble::as_tibble(
    Reduce(
      dplyr::bind_rows,
      lapply(
        mobster::Clusters(x) %>% dplyr::pull(cluster) %>% unique,
        dynamic
      )
    )
  )


}

# x ~ mobster fit
# cutoff_lv ~ named vector with the cutoffs to extract the latent variables of each cluster
# chromosomes ~ chromosomes that have been used for this analysis (could be a subset of the genome)
heuristic_mobster_clustering_enriched_for_chromsome_locations = function(x, cutoff_lv, chromosomes)
{
  stopifnot(all(cutoff_lv >= 0 & cutoff_lv <= 1))
  # stopifnot(N_min >= 0)
  stopifnot(inherits(x, 'dbpmm'))

  # Extract the location counts per cluster
  location_counter =
    lapply(
      mobster::Clusters(x) %>% dplyr::pull(cluster) %>% unique,
      function(cluster)
      {
        # Clusters
        assignments = mobster::Clusters(x, cutoff_assignment = cutoff_lv[cluster]) %>%
          dplyr::filter(cluster == !!cluster)

        tb = table(assignments$chr)

        # Count entries - explicit account for missing values (old school)
        missing_entries = setdiff(chromosomes, names(tb))
        entries = c(unlist(tb) %>% as.numeric(), rep(0, length(missing_entries)))
        names_entries = c(names(tb), missing_entries)

        df = data.frame(
          matrix(entries, ncol = length(entries), nrow = 1),
          stringsAsFactors = F
        )
        colnames(df) = names_entries

        df$cluster = cluster

        df
      })

  # These are the counts of the obersved mutations' chromosomes
  location_counter = Reduce(dplyr::bind_rows, location_counter) %>%
    as_tibble

  # Inspect each cluster
  evaluation =
    lapply(
      1:nrow(location_counter),
      function(i)
      {
        cluster = location_counter$cluster[i]

        # Extract cuonts
        mapping = location_counter %>%
          dplyr::select(-cluster) %>%
          dplyr::filter(row_number() == i)

        mapping = sort(mapping %>% as_vector(), decreasing = TRUE)

        # 1% of the total counts ~ what is less than that is not counted
        # but we do not above imposing > 10 mutations (10 is the minimum evidence
        # we decide we believe and that's it)
        min_cutoff_mapping = sum(mapping) * 0.01
        if(min_cutoff_mapping > 10) min_cutoff_mapping = 10

        which_above = which(mapping >= min_cutoff_mapping)

        # Special case, no CNA to match, we PASS by default
        if(length(which_above) == 0) {
          return(
            data.frame(
              cluster = cluster,
              prob = NA,
              p_20p = NA,
              CLUSTERED = FALSE,
              stringsAsFactors = F
            )
          )
        }

        mapping = mapping[which_above]

        # The mass are the total counts and the cumulative ones
        total_mass = sum(mapping)
        cumulative_mass = cumsum(mapping %>% as.numeric())

        # Special case: one chromosome is default enrichment
        if(length(cumulative_mass) == 1) {
          return(
            data.frame(
              cluster = cluster,
              prob = 1,
              p_20p = 0,
              CLUSTERED = TRUE,
              stringsAsFactors = F
            )
          )
        }

        # We compute now 60% of the total mass, and how many
        # cumulative counts are above that value, selecting the first
        # one that flags true
        index_prob = sum(cumulative_mass < total_mass * 0.6) + 1

        # We compare it to 20% of the required chromsome counts
        # If index_prob is above index_20p, then >60% of the observations
        # are occupying less than 20$ of the genome analysed, so they
        # are worth checking if they are a CNA-related event. These
        # cutoffs are arbitrary and chosen based on testing on WGS data.
        index_20p = length(which_above) * .2

        data.frame(
          cluster = cluster,
          prob = index_prob,
          p_20p = index_20p,
          CLUSTERED = index_prob > index_20p,
          stringsAsFactors = F
        )
      })

  # These are the tests for genome locations
  location_counter_table = Reduce(dplyr::bind_rows, evaluation) %>%
    as_tibble

  return(location_counter %>%
           dplyr::full_join(location_counter_table, by = 'cluster')
  )
}

# determines the sample purity without CNA data
heuristic_mobster_purity_estimation_nocna = function(x)
{
  L = sapply(
    mobster::Clusters(x) %>% dplyr::pull(cluster) %>% unique,
    function(cluster)
    {
      # We just assume everything is diploid, therefore 2 * the Beta peak of the clonal cluster
      estimated_tumour_purity = x$Clusters %>%
        filter(cluster %in% !!cluster, type == 'Mean') %>%
        pull(fit.value) %>% mean * 2
    }
  )

  # Some values might be  >1
  data.frame(
    cluster = mobster::Clusters(x) %>% dplyr::pull(cluster) %>% unique,
    purity = L,
    stringsAsFactors = FALSE
    )
}



