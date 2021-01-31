# Downloaded example PCAWG
x = readRDS("~/Downloads/evoverse_dataqc.rds")

require(tidyverse)

snvs = x$cnaqc$snvs %>% select(chr, from, to, ref, alt, NV, DP, VAF, gene, is_driver, driver_label)
cna = x$cnaqc$cna %>% select(chr, from, to, Major, minor)
purity = x$cnaqc$purity
reference = x$cnaqc$reference_genome

require(evoverse)

options(easypar.parallel = FALSE)

data_qc = pipeline_qc_copynumbercalls(
  mutations = snvs,
  cna = cna,
  purity = purity,
  reference = 'hg19',
  description = "Data QC example",
  smooth = TRUE,
  matching_epsilon_peaks = 0.035,
  ccf_method = 'ROUGH',
  peak_method = 'closest',
  min_CCF = 0.1,
  only_SNVs = TRUE
)

saveRDS(data_qc, file = "R/nobuild/dataqc.rds")

ggsave("R/nobuild/dataqc.pdf", plot = plot(data_qc), width = 10, height = 12)

# Timing deconvolution

timing_fit = pipeline_chromosome_timing(
  data_qc,
  karyotypes = c('2:0', '2:1', '2:2'),
  min_muts = 50,
  description = "Chromosomal timing sample",
  N_max = 15000,
  enforce_QC_PASS = FALSE,
  auto_setup = "FAST"
)

saveRDS(timing_fit, file = "R/nobuild/timing_fit.rds")

ggsave("R/nobuild/timingfitc.pdf", plot = plot(timing_fit), width = 10, height = 10)

# VAF deconvolution

vaf_fit = pipeline_subclonal_deconvolution_VAF_karyotype(
  data_qc,
  min_muts = 50,
  description = "VAF deconvolution by karyotype",
  N_max = 15000,
  enforce_QC_PASS = FALSE,
  BMix_entropy = TRUE,
  auto_setup = "FAST"
)

saveRDS(vaf_fit, file = "R/nobuild/vaf_fit.rds")

ggsave("R/nobuild/vaf_fit.pdf", plot = plot(vaf_fit), width = 10, height = 10)

# CCF deconvolution

ccf_fit = pipeline_subclonal_deconvolution_CCF(
  data_qc,
  karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
  BMix_entropy = TRUE,
  enforce_QC_PASS = FALSE,
  min_muts = 50,
  description = "CCF deconvolution",
  N_max = 15000,
  auto_setup = "FAST"
)

saveRDS(ccf_fit, file = "R/nobuild/ccf_fit.rds")

ggsave("R/nobuild/ccf_fit.pdf", plot = plot(ccf_fit), width = 10, height = 10)
