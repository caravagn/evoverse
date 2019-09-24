load('~/Downloads/Archive/data/Mutect_Set07_Diploid.RData', verbose = T)

mutations = readxl::read_excel('~/Downloads/Archive/SM_tables/Supplementary_Table_S1_MSeq_calls.xlsx', sheet = 2)

cna = readxl::read_excel('~/Downloads/Archive/SM_tables/Supplementary_Table_S1_MSeq_calls.xlsx', sheet = 3) 
colnames(cna)[c(2,4)] = c('from', 'to')
colnames(cna) = gsub('_Minor', '\\.minor', colnames(cna))
colnames(cna) = gsub('_Major', '\\.Major', colnames(cna))
cna$chr = paste0('chr', cna$chr)

clean_mutations = mutations %>%
  select(-ends_with('projected'), -ends_with('cluster'),-starts_with('MOBSTER')) %>%
  mutate(from = as.numeric(from), to = as.numeric(to))
colnames(clean_mutations) = gsub('_original', '', colnames(clean_mutations))

mutations_list = cna_list = NULL

for(s in dataset$samples)
{
  nw = clean_mutations %>% select(chr,
                           from,
                           to,
                           ref,
                           alt,
                           gene,
                           starts_with(!!s),
                           cosmic,
                           `function.`,
                           mutlocation,
                           region)

  colnames(nw) = gsub(paste0(s, '\\.'), '', colnames(nw))
  colnames(nw)[7] = 'vaf_normal'
  
  mutations_list = append(mutations_list, list(nw))

  # 
  nw = cna %>% select(chr,
                      from,
                      to,
                      nloci,
                      starts_with(!!s))
  
  colnames(nw) = gsub(paste0(s, '\\.'), '', colnames(nw))

  cna_list = append(cna_list, list(nw))
}

names(cna_list) = names(mutations_list) = dataset$samples

save(clean_mutations, file = '../Set7.RData')

require(mobster)
require(mvmobster)

x = mvmobster::dataset(
  mutations = mutations_list,
  segments = cna_list,
  samples = dataset$samples,
  purity = dataset$purity
)

library(CNAqc)

f = plot(x)
f[[3]]

require(ggpubr)
fg = ggarrange(plotlist = f, ncol = 2, nrow = 3)

plot_sample(x, sample = 'Set7_55')

qc = CNAqc::analyze_peaks(x$CNAqc$Set7_62)
plot_peaks_analysis(qc)

qc$peaks_analysis$plots$`1:1` + ylim(0, 5)

y = analyze_mobster(x, parallel = F, epsilon = 1e-5, K = c(1, 2))
save(y, file = '../y_Set7.RData')


f = plot.mbst_data(y, clusters = 'MOBSTER')
f[[3]]


