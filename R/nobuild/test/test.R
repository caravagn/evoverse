x = evoverse::pipeline_qc_copynumbercalls(
  mutations = CNAqc::example_dataset_CNAqc$snvs,
  cna = CNAqc::example_dataset_CNAqc$cna,
  purity = CNAqc::example_dataset_CNAqc$purity
)

ccf = evoverse::pipeline_chromosome_timing(x, auto_setup = 'FAST')

vaf = evoverse::pipeline_chromosome_timing(x, auto_setup = 'FAST')

ccf = pipeline_subclonal_deconvolution_CCF(x)


require(devtools)

install("../CNAqc/");load_all("~/Documents/GitHub/CNAqc")

load_all(".")

load('/Volumes/Data/Dropbox/160711_HMFregCPCT_FR12244551_FR12244640_CPCT02020325_cnaqc_object.RData')

x = evoverse::pipeline_qc_copynumbercalls(
  mutations = x$snvs,
  cna = x$cna %>% select(-n, -segment_id),
  purity = x$purity,
  only_SNVs = T,
  min_CCF = 0.15
)

plot_data_histogram(x$cnaqc)

xvaf = evoverse::pipeline_subclonal_deconvolution_VAF_karyotype(x, auto_setup = "FAST")
xccf = evoverse::pipeline_subclonal_deconvolution_CCF(x, auto_setup = "FAST")

pdf("a.pdf", width = 10, height = 9)
plot(xvaf)
plot(xccf)
dev.off()


Reduce(bind_rows, lapply(xvaf$mobster, mobster::Clusters))


dnds = lapply(xvaf$mobster, function(x) {
  dnds_stats = dnds(mobster::Clusters(x),
                    gene_list = mobster::cancer_genes_dnds$Tarabichi_drivers)
})


dnds_stats = dnds(mobster::Clusters(xvaf$mobster$`2:1`))
print(dnds_stats$plot)

dnds_multi = dnds(
  Reduce(bind_rows,
         lapply(xvaf$mobster %>% names, function(x)
           mobster::Clusters(xvaf$mobster[[x]]) %>%
             select(chr, from, ref, alt, cluster) %>%
             mutate(sample = x)
           )
         ),
  mapping = c(`C3` = 'Non-tail',
              `C2` = 'Non-tail',
              `C1` = 'Non-tail',  # Pool together all clonal mutations
              `Tail` = 'Tail'     # Pool together all tail mutations),
  )
)

print(dnds_multi$plot)


install.packages('ActiveDriverWGS', dependencies = T)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(ActiveDriverWGS)

data("cll_mutations")
head(cll_mutations)
data("cancer_genes")
head(cancer_genes)
some_genes = c("ATM", "MYD88", "NOTCH1")

results = ActiveDriverWGS(mutations = cll_mutations,
                          elements = cancer_genes[cancer_genes$id %in% some_genes,],
                          sites = cancer_gene_sites)

all_mutations = Reduce(bind_rows,
       lapply(xvaf$mobster %>% names, function(x)
         mobster::Clusters(xvaf$mobster[[x]]) %>%
           select(chr, from, to, ref, alt, cluster, gene) %>%
           mutate(sample = x, pos1 = from, pos2 = to, patient =x)
       )
)

some_genes = mobster::cancer_genes_dnds$Martincorena_drivers

results = ActiveDriverWGS(mutations = all_mutations,
                          # elements = cancer_genes[cancer_genes$id %in% some_genes,],
                          # sites = cancer_gene_sites,
                          ref_genome = 'hg38')

