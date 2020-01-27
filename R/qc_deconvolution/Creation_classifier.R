##############################################################################################################################
# Get new versions from private repos
devtools::install_github('caravagn/mobster', auth_token = '0df3389381cd57ef4636cb5aa8e9f81f78d14e81', ref = 'development')
library(tidyverse)

##############################################################################################################################
# DAVROS: in folder fits, extracts a table of the model parameters
extract_covariates = function(p)
{
  cat("\n", p)

  if(file.exists(paste0(p, "/table_deconvolution.csv"))) {
    cat('CACHE')
    return(2)
  }

  load(paste0(p, '/HW_deconvolution.RData'))


  Reduce(
    bind_rows,
    lapply(
    names(karyotypes_analysis),
    function(k)
    {
      if(all(is.null(karyotypes_analysis[[k]]$MOBSTER_fit))) return(NULL)

      mobster::to_string(karyotypes_analysis[[k]]$MOBSTER_fit$best) %>%
        mutate(
          sample = p,
          karyotype = k
        )
    })
  ) %>%
    write_tsv(paste0(p, "/table_deconvolution.csv"))

  cat('DONE ', p)

}

# extract_covariates('160711_HMFregCPCT_FR12244551_FR12244640_CPCT02020325')


require(tidyverse)

lf = list.dirs('.')
lf = lf[lf != '']

easypar::run(FUN = extract_covariates, PARAMS = lapply(lf, list), parallel = F)


##############################################################################################################################
# DAVROS: Assemble the results from the above code in one table
samples = list.dirs('.')
samples = samples[samples != '']
library(tidyverse)

rd = function(p) {  read_tsv(paste0(p, '/table_deconvolution.csv'))  }
L = easypar::run(FUN=rd, PARAMS=lapply(samples, list), parallel = F)
L = Reduce(bind_rows, L)

write_tsv(L, "Parameters_fits.tsv")
#



##############################################################################################################################
# Locally -- multinomial classifier

# library(ggforce)
# library(ggfortify)

parameters_models =  read_tsv("Parameters_fits.tsv")
parameters_models$sample = gsub('\\./', '', parameters_models$sample)

annotated_fits =  read_tsv("training_set.tsv") %>%
  mutate(
    karyotype = paste(karyotype),
    karyotype = case_when(
      karyotype == "01:00:00" ~ '1:0',
      karyotype == "01:01:00" ~ '1:1',
      karyotype == "02:01:00" ~ '2:1',
      karyotype == "02:00:00" ~ '2:0',
      karyotype == "02:02:00" ~ '2:2'
      )
  ) %>%
  rename(Training_Set_QC = QC)

# Covariates to predict
parameters_models =
  parameters_models %>%
  select(
    starts_with("Mean_C"),
    starts_with("Variance_C"),
    N,
    starts_with("N_C"),
    starts_with("pi_"),
    Shape_Tail,
    entropy,
    reduced.entropy,
    sse_total,
    sse_0_0.1,
    sse_0.1_0.2,
    sse_0.6_0.7,
    sse_0.7_0.8,
    sse_0.8_0.9,
    sse_0.9_1,
    sample,
    karyotype
  )  %>%
  full_join(annotated_fits) %>%
  filter(!is.na(Training_Set_QC))


parameters_models %>%
  ggplot(aes(sse_0_0.1, fill = Training_Set_QC)) +
  geom_histogram()

# load the package
library(VGAM)

df_params = parameters_models %>%
  select(starts_with('sse'), starts_with("Variance_C"),
         Shape_Tail,
         entropy,
         reduced.entropy, Training_Set_QC, sample, karyotype) %>%
  filter(!is.na(sse_total))

CASES = table(df_params$Training_Set_QC)
SPLITCASES = round(CASES * .7)

# f  p r -- > random preserved proportions
df_training =
  bind_rows(
    df_params %>% filter(Training_Set_QC == 'f') %>% sample_n(SPLITCASES['f']),
    df_params %>% filter(Training_Set_QC == 'p') %>% sample_n(SPLITCASES['p']),
    df_params %>% filter(Training_Set_QC == 'r') %>% sample_n(SPLITCASES['r'])
  ) %>%
  select(sample, karyotype) %>%
  mutate(VGLM_TRAINING = TRUE)

df_params = full_join(df_params, df_training) %>%
  mutate(VGLM_TRAINING = ifelse(is.na(VGLM_TRAINING), FALSE, VGLM_TRAINING))

training_partition = df_params %>%
  filter(VGLM_TRAINING) %>%
  select(-karyotype, -sample, -VGLM_TRAINING)

test_partition = df_params %>%
  filter(!VGLM_TRAINING) %>%
  select(-karyotype, -sample, -VGLM_TRAINING)

# fit model
fit <- VGAM::vglm(Training_Set_QC~.,
                  family = VGAM::multinomial,
                  data = training_partition)

# summarize the fit
summary(fit)

# make predictions
# predicitons = NULL
#
# for(j in 1:nrow(test_partition))
# {
#   class_prob = VGAM::predict(
#     fit,
#     data_test,
#     type = "response")
#
#   print(j)
#   print(class_prob)
#
#   predictions = rbind(predictions, as_tibble(class_prob))
# }

probabilities = VGAM::predict(
  fit,
  test_partition,
  type = "response")

test_partition$prediction = colnames(probabilities)[apply(probabilities, 1, which.max)]
sum(test_partition$Training_Set_QC == test_partition$prediction)/nrow(test_partition)

# summarize accuracy


# Logistic classifier p vs NONp
parameters_models =  read_tsv("Parameters_fits.tsv")
parameters_models$sample = gsub('\\./', '', parameters_models$sample)

annotated_fits =  read_tsv("training_set.tsv") %>%
  mutate(
    karyotype = paste(karyotype),
    karyotype = case_when(
      karyotype == "01:00:00" ~ '1:0',
      karyotype == "01:01:00" ~ '1:1',
      karyotype == "02:01:00" ~ '2:1',
      karyotype == "02:00:00" ~ '2:0',
      karyotype == "02:02:00" ~ '2:2'
    )
  ) %>%
  rename(Training_Set_QC = QC)

# Covariates to predict
parameters_models =
  parameters_models %>%
  select(
    starts_with("Mean_C"),
    starts_with("Variance_C"),
    N,
    starts_with("N_C"),
    starts_with("pi_"),
    Shape_Tail,
    entropy,
    reduced.entropy,
    sse_total,
    sse_0_0.1,
    sse_0.1_0.2,
    sse_0.6_0.7,
    sse_0.7_0.8,
    sse_0.8_0.9,
    sse_0.9_1,
    sample,
    karyotype
  )  %>%
  full_join(annotated_fits) %>%
  filter(!is.na(Training_Set_QC))


parameters_models %>%
  ggplot(aes(x = sse_total, y = sse_0_0.1, color = Training_Set_QC)) +
  geom_jitter() +
  facet_wrap(~Training_Set_QC)

parameters_models %>%
  ggplot(aes(x = sse_0.1_0.2, y = sse_0_0.1, color = Training_Set_QC)) +
  geom_jitter() +
  facet_wrap(~Training_Set_QC)



load('~/a.RData')
M

genes_list = mobster::cancer_genes_dnds$Martincorena_drivers

mobster:::extract_dnds_input(x=M, genes_list = genes_list) %>% pull(gene)

data('fit_example')
x = fit_example$best


x_domain = seq(0, 1, 0.001)
domain = mobster::Clusters_denovo(x, data.frame(VAF = x_domain))

densities = mobster:::template_density(
  x,
  x.axis = seq(1e-3, 1-1e-3, by = 1e-3), # Restricted for numerical errors
  binwidth = 1e-3,
  reduce = TRUE) %>%
  spread(cluster, y)

