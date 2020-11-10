x = evoverse::pipeline_qc_copynumbercalls(
  mutations = CNAqc::example_dataset_CNAqc$snvs,
  cna = CNAqc::example_dataset_CNAqc$cna,
  purity = CNAqc::example_dataset_CNAqc$purity
)

ccf = evoverse::pipeline_chromosome_timing(x, auto_setup = 'FAST')
vaf = evoverse::pipeline_chromosome_timing(x, auto_setup = 'FAST')

ccf = pipeline_subclonal_deconvolution_CCF(x)


