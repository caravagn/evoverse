data('example_mvmobster')
require(tidyverse)

m1 = example_mvmobster$muts %>% select(region, gene, chr, from, to, ref, alt, Set6_42.VAF, Set6_42.DP, Set6_42.NV, patient) %>% filter(Set6_42.DP > 0)
m2 = example_mvmobster$muts %>% select(region, gene, chr, from, to, ref, alt, Set6_44.VAF, Set6_44.DP, Set6_44.NV, patient) %>% filter(Set6_44.DP > 0)

colnames(m1) = c('region', 'gene', 'chr', 'from', 'to', 'ref', 'alt', 'VAF', 'DP', 'NV', 'patient')
colnames(m2) = c('region', 'gene', 'chr', 'from', 'to', 'ref', 'alt', 'VAF', 'DP', 'NV', 'patient')

s1 = example_mvmobster$segments %>% select(chr, from, to, Set6_42.minor, Set6_42.Major)
s2 = example_mvmobster$segments %>% select(chr, from, to, Set6_44.minor, Set6_44.Major)

colnames(s1) = c('chr', 'from', 'to', 'minor', 'Major')
colnames(s2) = c('chr', 'from', 'to', 'minor', 'Major')

purity = c(.8, .8)

mutations = list(`Set42` = m1, `Set44` = m2)
segments = list(`Set42` = s1, `Set44` = s2)
names(purity) = c('Set42', 'Set44')
samples = names(purity)

