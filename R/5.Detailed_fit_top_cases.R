library(tidyverse)

results = readr::read_csv(file = 'PCWAG_Summary_table.csv')

# Load table
PCWAG = read_tsv('consensus.20170217.purity.ploidy.txt')

results = results %>%
  left_join(PCWAG, by = 'samplename')

results = results %>% mutate(tail = ifelse(tail, "With tail", "Without tail"))
results$tail = factor(results$tail, levels = c("With tail", "Without tail"))

results = results %>%
  mutate(
    subclone = ifelse(
      K_beta > 1 &
        pi_C2 <= pi_C1 &
        Mean_C1 < 0.5,
      TRUE,
      FALSE
    ))



timon_plots = list.files('example_comparative/')
timon_plots = gsub(timon_plots, pattern = '_VAF_histrogram_all.pdf', replacement = '')
timon_plots = gsub(timon_plots, pattern = '_VAF_histrogram.pdf', replacement = '')
timon_plots = timon_plots %>% unique

with_timon_plots = intersect(results$samplename, timon_plots)

cases = results %>% mutate(
  rank = purity * coverage
) %>%
  select(-1, -2, -starts_with('Variance')) %>%
  arrange(desc(rank), N_Tail) %>%
  filter(samplename %in% with_timon_plots) %>%
  filter(row_number() < 20)

zz = sapply(cases$samplename %>% unique, FUN = copy_from_davros, DATA = TRUE, FOLDER = '~/Documents/Andrea/')
zz = de_novo_plots(cases)

assemble_image(cases, "Top 20 quality rank - with Timon plots", 'Top_20_with_Timon_plot.pdf', add_timon = TRUE)

