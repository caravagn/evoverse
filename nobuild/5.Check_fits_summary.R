all_files = read_tsv('consensus.20170217.purity.ploidy.txt') %>%
  mutate(
    tsv_file =
      file.exists(
        paste0(
          'fits/', samplename, '-injected_CN_ICL.tsv'
        )
      )
  )

table(all_files$tsv_file)

