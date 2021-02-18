require(tidyverse)

# Example data in the CNAqc package - input formats are discussed
# here: https://caravagnalab.github.io/CNAqc/articles/a1_Introduction.html
snvs = CNAqc::example_dataset_CNAqc$snvs
cna = CNAqc::example_dataset_CNAqc$cna
purity = CNAqc::example_dataset_CNAqc$purity
reference = CNAqc::example_dataset_CNAqc$reference

require(evoverse)

options(easypar.parallel = FALSE)

data_qc = pipeline_qc_copynumbercalls(
  mutations = snvs,
  cna = cna,
  purity = purity,
  reference = 'GRCh38',
  description = "Data QC example",
  smooth = TRUE,
  matching_epsilon_peaks = 0.035,
  ccf_method = 'ENTROPY',
  peak_method = 'closest',
  min_CCF = 0.1,
  only_SNVs = TRUE
)

saveRDS(data_qc, file = "dataqc.rds")

ggsave("dataqc.pdf", plot = plot(data_qc), width = 10, height = 12)

# VAF deconvolution

vaf_fit = pipeline_subclonal_deconvolution_VAF_karyotype(
  data_qc,
  min_muts = 50,
  karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
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
