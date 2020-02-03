load("~/Documents/Davros/Hartwig_analysis/fits/190913_HMFregCPCT_FR16672425_FR13923148_CPCT02020940/HW_input_data.RData")
hw_data

load("~/Downloads//HW_input_data.RData")

# CNAqc::plot_icon_CNA(hw_data)
# CNAqc::plot_segments(hw_data, circular = T)
#
# hw_data = smooth_segments(hw_data, maximum_distance = 5e6)
#
# hw_data = detect_wg_overfragmentation(hw_data)
# hw_data = detect_arm_overfragmentation(hw_data)
#
# plot_arm_fragmentation(hw_data)
# plot_segment_size_distribution(hw_data)
# plot_segments(hw_data)
# plot_smoothing(hw_data)


ct_deconvolution = evoverse::pipeline_chromosome_timing(
  mutations = hw_data$snvs,
  cna = hw_data$cna,
  purity = hw_data$purity,
  N_max = 3000,
  description = "My amazing dataset!",
  auto_setup = 'FAST'
)

sdv_deconvolution = evoverse::pipeline_subclonal_deconvolution_VAF(
  mutations = hw_data$snvs,
  cna = hw_data$cna,
  purity = hw_data$purity,
  N_max = 3000,
  description = "160704_HMFregCPCT_FR12244557_FR12244595_CPCT02110002",
  auto_setup = 'FAST'
)


sdc_deconvolution = evoverse::pipeline_subclonal_deconvolution_CCF(
  mutations = hw_data$snvs,
  cna = hw_data$cna,
  purity = hw_data$purity,
  N_max = 3000,
  description = "160704_HMFregCPCT_FR12244557_FR12244595_CPCT02110002",
  auto_setup = 'FAST'
)

lf = list.dirs("~/Documents/Davros/Hartwig_analysis/fits", full.names = F)
lf = lf[lf != '']

expand.grid(`karyotype` = c('1:0', '1:1', '2:1', '2:0', '2:2'), `sample` = lf) %>%
  as_tibble() %>%
  write_tsv('~/Documents/Github/evoverse/R/qc_deconvolution/training_set.tsv')

