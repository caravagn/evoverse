x = readr::read_tsv('~/Downloads/test.vaf', col_names = F)
cols = readr::read_tsv(
  '~/Downloads/final_consensus_passonly.snv_mnv_indel.icgc.public.maf',
  col_names = F,
  n_max = 1
)
colnames(x) = cols[1, ]

ggplot(x, aes(as.numeric(i_VAF))) + geom_histogram(binwidth = 0.01)

x = x %>% select(
  Chromosome,
  Start_position,
  End_position,
  Reference_Allele,
  Tumor_Seq_Allele1,
  t_alt_count,
  t_ref_count
)
colnames(x) = c('chr', 'from', 'to', 'ref', 'alt', 'NV', 'NR')
x$DP = x$NV + x$NR
x$VAF = x$NV / x$DP

x$chr = paste0('chr', x$chr)


segm = readr::read_tsv(
  '~/Downloads/0b29c893-03bf-4131-b192-c14a2788d411.consensus.20170119.somatic.cna.annotated.txt'
) %>%
  select(chromosome, start, end, major_cn, minor_cn, battenberg_frac1_A) %>%
  filter(battenberg_frac1_A == 1)
colnames(segm) = c('chr', 'from', 'to', 'Major', 'minor', 'CCF')

cnqx = CNAqc::init(x, segm, '.99', ref = 'hg19') %>% CNAqc::subset_snvs()
CNAqc::plot_data_histogram(cnqx)


PCAWG_mutations = readr::read_tsv("~/Downloads/final_consensus_passonly.snv_mnv_indel.icgc.public.maf")

setwd('Research Projects/2020. CNAqc/pcawg/good/')
files = list.files()
files = files[-11]
files = gsub('.png', '', files)

for (s in files)
{
  sm = PCAWG_mutations %>% filter(Tumor_Sample_Barcode == s) %>% 
    select(Chromosome, ends_with('position'), ends_with('Allele'),Tumor_Seq_Allele1, t_alt_count, t_ref_count)
  
  colnames(sm) = c('chr', 'from', 'to', 'ref', 'alt', 'NV','NR')
  sm$DP=sm$NR+sm$NV
  sm$VAF = sm$NV/sm$DP
  sm$chr=paste0('chr',sm$chr)
  
  segm = readr::read_tsv(
    paste0('~/Downloads/', s, '.consensus.20170119.somatic.cna.annotated.txt')
  ) %>%
    select(chromosome,
           start,
           end,
           major_cn,
           minor_cn,
           battenberg_frac1_A) %>%
    filter(battenberg_frac1_A == 1)
  colnames(segm) = c('chr', 'from', 'to', 'Major', 'minor', 'CCF')
  
  if(nrow(segm) == 0) next
  
  
  cnqx = CNAqc::init(sm, segm, '.99', ref = 'hg19') %>% CNAqc::subset_snvs()
  
  CNAqc::plot_data_histogram(cnqx)+labs(title=s) +
    ggsave(paste0(s,'_n.png'))
  
}
