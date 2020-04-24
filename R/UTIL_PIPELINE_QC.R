# Get a QC table
get_pqc_cna_QC_table = function(x)
{
  all_karyptypes = x$peaks_analysis$fits %>% names

  if(all(is.null(x$peaks_analysis)))
  {
    peaks_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'Peaks',
      stringsAsFactors = FALSE
    )
  }
  else
  {
    peaks_QC = x$peaks_analysis$matches %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
      dplyr::distinct(karyotype, QC, .keep_all = T) %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(
        value = 1,
        lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        label = karyotype,
        type = 'Peaks')
  }

  if(all(is.null(x$CCF_estimates)))
  {
    CCF_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'CCF',
      stringsAsFactors = FALSE
    )
  }
  else
  {
    CCF_QC = Reduce(dplyr::bind_rows, lapply(x$CCF_estimates, function(x) x$QC_table)) %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(
        value = 1,
        lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        label = karyotype,
        type = 'CCF')
  }

  QC_table = dplyr::bind_rows(peaks_QC, CCF_QC)
  QC_table$karyotype = factor(QC_table$karyotype, all_karyptypes)
  QC_table$type = factor(QC_table$type, levels = c('Peaks', 'CCF'))

  QC_table %>%
    dplyr::mutate(
      QC = ifelse(!is.na(QC) & QC == "NA", NA, QC)
    )
}

