setwd("~/Documents/Github/mvMOBSTER")
R.utils::sourceDirectory("old", modifiedOnly=FALSE)
R.utils::sourceDirectory("R", modifiedOnly=FALSE)

data('example_mvmobster')
require(tidyverse)
require(pio)

m1 = example_mvmobster$muts %>% select(chr, from, to, ref, alt, region, gene, Set6_42.VAF, Set6_42.DP, Set6_42.NV, patient) %>% filter(Set6_42.DP > 0)
m2 = example_mvmobster$muts %>% select(chr, from, to, ref, alt, region, gene, Set6_44.VAF, Set6_44.DP, Set6_44.NV, patient) %>% filter(Set6_44.DP > 0)

colnames(m1) = c('chr', 'from', 'to', 'ref', 'alt', 'region', 'gene', 'VAF', 'DP', 'NV', 'patient')
colnames(m2) = c('chr', 'from', 'to', 'ref', 'alt', 'region', 'gene', 'VAF', 'DP', 'NV', 'patient')

s1 = example_mvmobster$segments %>% select(chr, from, to, Set6_42.minor, Set6_42.Major) %>% as_tibble()
s2 = example_mvmobster$segments %>% select(chr, from, to, Set6_44.minor, Set6_44.Major) %>% as_tibble()

colnames(s1) = c('chr', 'from', 'to', 'minor', 'Major')
colnames(s2) = c('chr', 'from', 'to', 'minor', 'Major')

purity = c(.8, .8)

mutations = list(`Set42` = m1, `Set44` = m2)
segments = list(`Set42` = s1, `Set44` = s2)
names(purity) = c('Set42', 'Set44')
samples = names(purity)

x = dataset(
  mutations,
  segments,
  samples,
  purity
)

VAF(x)
plot(x)

x = analyze_mobster(x, parallel = F, K = 1,  init = 'random')

plot(x, clusters = 'MOBSTER')

x_notail = filter_tails(x)
x_notail = analyze_VIBER(x_notail)

plot(x_notail, clusters = 'VIBER')

x = analyze_VIBER(x)
plot(x, clusters = 'VIBER')


#
load('../test.dbpmm/Real Data/[Used] Koerber_et_al/data/H043-4PGF.RData')
cp = dataset$primary_CNA %>%
  select(chr, from, to, genotype) %>%
  filter(!is.na(genotype)) %>%
  separate(col = 'genotype', sep = ':', into = c('minor', 'Major')) %>%
  filter(minor != 'sub', Major != 'sub') %>%
  as_tibble()

cr = dataset$relapse_CNA %>%
  select(chr, from, to, genotype) %>%
  filter(!is.na(genotype)) %>%
  separate(col = 'genotype', sep = ':', into = c('minor', 'Major')) %>%
  filter(minor != 'sub', Major != 'sub') %>%
  as_tibble()


example_inputs_mutations_wide = example_mvmobster$muts %>% select(chr, from, to, ref, alt, region, gene,
                                  Set6_42.VAF, Set6_42.DP, Set6_42.NV,
                                  Set6_44.VAF, Set6_44.DP, Set6_44.NV,
                                  patient) %>% as_tibble()

example_inputs_segments_wide = example_mvmobster$segments %>% select(chr, from, to,
                                                                     Set6_42.minor, Set6_42.Major,
                                                                     Set6_44.minor, Set6_44.Major) %>% as_tibble()

example_input_formats = list(
  wide_mutations = example_inputs_mutations_wide,
  wide_segments = example_inputs_segments_wide,
  long_mutations = mutations,
  long_segments = segments
)

usethis::use_data(example_input_formats, overwrite = T)

example_inputs_mutations_wide = example_mvmobster$muts %>% select(chr, from, to, ref, alt, region, gene,
                                                                  Set6_42.VAF, Set6_42.DP, Set6_42.NV,
                                                                  Set6_44.VAF, Set6_44.DP, Set6_44.NV,
                                                                  patient)


