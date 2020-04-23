library(ggforce)
library(ggfortify)
library(tidyverse)



dfQC <- read.csv("annotated_training_set.csv", stringsAsFactors = F)

dfQCf <- dfQC %>%
  filter(N_C1 > 2, N_C2 > 2) %>%
  group_by(id, karyotype) %>%
  mutate(ratioM = abs((Mean_C1 / Mean_C2) - 2),
         ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
         ratioN = (max(N_C1, N_C2) +1) / (min(N_C2, N_C1) + 1),
         minpi = min(pi_C1, pi_C2)) %>%
  ungroup() %>%
  filter(PF != "")


dfpca <- dfQCf %>%
  select(N, ratioM, ratioV, ratioN, reduced.entropy, tail)
pca <- prcomp(dfpca, scale = TRUE)
dfpca <- bind_cols(dfpca, as.data.frame(pca$x))

(g1 <- autoplot(pca, data = dfQCf, alpha = 0.5,
                scale = 0, loadings = TRUE, loadings.label = T,
                x = 1, y = 2, colour = "PF"))
(g1 <- autoplot(pca, data = dfQCf, alpha = 0.5,
                scale = 0,
                x = 1, y = 2, colour = "PF"))

(g <- dfQCf %>%
    select(N, ratioM, ratioV, Variance_C1, Variance_C2,
           ratioN, reduced.entropy, PF, tail, minpi) %>%
    gather(key = "param", value = "value", -PF) %>%
    ggplot(aes(x = value, fill = PF)) +
    geom_histogram(bins = 50) +
    facet_wrap(~param, scale = "free"))
cowplot::save_plot("histogram.pdf", g, base_height = 10, base_width = 10)

(g <- dfQCf %>%
    select(N, ratioM, ratioV, Variance_C1, Variance_C2,
           ratioN, reduced.entropy, PF, tail, minpi) %>%
    gather(key = "param", value = "value", -PF) %>%
    ggplot(aes(x = value, fill = PF)) +
    geom_histogram(bins = 50) +
    scale_x_log10() +
    facet_wrap(~param, scale = "free"))
cowplot::save_plot("histogramlog.pdf", g, base_height = 10, base_width = 10)




set.seed(123)

dfQCf <- dfQC %>%
  filter(N_C1 > 2, N_C2 > 2) %>%
  group_by(id, karyotype) %>%
  mutate(ratioM = abs((Mean_C1 / Mean_C2) - 2),
         ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
         ratioN = (max(N_C1, N_C2) +1) / (min(N_C2, N_C1) + 1),
         minpi = min(pi_C1, pi_C2)) %>%
  ungroup() %>%
  filter(PF != "")

dat <- dfQCf %>%
  select(minpi, ratioM, ratioV, ratioN,
         reduced.entropy, tail, PF, entropy) %>%
  mutate(PF = ifelse(PF == "P", 1, 0))
training.samples <- dat$PF %>%
  caret::createDataPartition(p = 0.8, list = FALSE)
train.data  <- dat[training.samples, ]
test.data <- dat[-training.samples, ]

model <- glm( PF ~., data = train.data, family = binomial)
saveRDS(object = model, file = "logsitic_model_timing.RData")
# Summarize the model
summary(model)
# Make predictions
probabilities <- model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Model accuracy
mean(predicted.classes == test.data$PF)


fulldat <- read_csv("timing_QC_table.csv") %>%
  #filter(N_C1 > 2, N_C2 > 2) %>%
  group_by(id, karyotype) %>%
  mutate(ratioM = abs((Mean_C1 / Mean_C2) - 2),
         ratioV = max(Variance_C1, Variance_C2) / min(Variance_C1, Variance_C2),
         ratioN = (max(N_C1, N_C2) +1) / (min(N_C2, N_C1) + 1),
         minpi = min(pi_C1, pi_C2)) %>%
  ungroup()

fulldatasetprob <- model %>% predict(fulldat, type = "response")
predicted.classes.fulldataset <- ifelse(fulldatasetprob > 0.5, "PASS", "FAIL")
table(predicted.classes.fulldataset)

fulldat$QC <- predicted.classes.fulldataset
fulldat$prob <- fulldatasetprob

write_csv(fulldat, "QC/QCresults.csv")
(g <- fulldat %>%
    ggplot(aes(x = prob)) +
    geom_histogram(bins = 100) +
    xlab("Probability of pass"))
save_plot("QC/prob_pass.pdf", g)
(g <- fulldat %>%
    ggplot(aes(x = prob)) +
    geom_histogram(bins = 100) +
    xlab("Probability of pass") +
    facet_wrap(~karyotype, scales = "free"))
save_plot("QC/prob_pass_kt.pdf", g)

(g1 <- fulldat %>%
  group_by(id) %>%
  summarize(prob = max(prob)) %>%
  ungroup() %>%
  arrange(prob) %>%
  mutate(n =  100 * (n() - 1:n()) / n()) %>%
  ggplot(aes(x = prob, y = n)) +
  geom_line() +
  xlab("Probability threshold") +
  ylab("% of samples retained (1+ karyotypes)"))
(g2 <- fulldat %>%
    group_by(id) %>%
    summarize(prob = max(prob)) %>%
    ungroup() %>%
    arrange(prob) %>%
    mutate(n =  (n() - 1:n())) %>%
    ggplot(aes(x = prob, y = n)) +
    geom_line() +
    xlab("Probability threshold") +
    ylab("# of samples retained (1+ karyotypes)"))
g <- plot_grid(g1, g2)
save_plot("QC/prob_threshold.pdf", g, base_aspect_ratio = 2)

fulldat %>%
  ggplot(aes(x = ratioN, y = N, col = QC)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

# Replot
plotdir <- "~/Documents/chrtiming/QCplots/"
basedir <- "~/re_gecip/shared_allGeCIPs/R_dm/data/COLORECTAL_timing/"

load("Timing_Status.RData")

for (i in 1:length(status)){
  if (status[[i]]$Timing_plot == FALSE){
    message("No timing results")
    next
  }
  if (round(i / 25) == (i / 25)){
    message(paste0(i, "/", length(status), " samples finished"))
  }
  myid <- status[[i]]$participant_id
  fits <- list.files(paste0(basedir, myid), pattern = "fit_and_data")
  plotlist <- list()
  j <- 1
  if (length(fits) == 0) {
    next
  }
  for (f in fits){
    load(paste0(basedir, myid, "/", f))
    mykaryotype <- str_split(f, "_")[[1]][5]
    QC <- filter(fulldat, id == myid, karyotype == mykaryotype)
    PF <- QC$QC
    prob <- round(QC$prob, 3)
    ggtitle1 <- paste0("ID = ", myid, ", Karyotype = ", mykaryotype)
    print(ggtitle1)
    ggsubtitle <- paste0(PF, " (p (PASS) = ", prob, ")")
    if (PF == "FAIL"){
      titlecol = "firebrick"
    } else {
      titlecol = "darkolivegreen"
    }
    plotlist[[j]] <- plot(plot$fit$best) + ggtitle(ggtitle1, subtitle = ggsubtitle) +
      theme(plot.title = element_text(color = titlecol),
            plot.subtitle = element_text(color = titlecol))
    j <- j + 1
  }
  g <- cowplot::plot_grid(plotlist = plotlist, ncol = 3)
  save_plot(paste0(plotdir, "persample/", myid, ".pdf"), g, base_height = 4, base_width = 10)
}




# Generate script for refitting
refitting <- fulldat %>%
  filter(QC == "FAIL") %>%
  filter(ratioN < 50) %>%
  distinct(id, karyotype) %>%
  rename(participant_id = id)

dfMSS <- read_delim("Timing_MSS_QCPASS.tsv", delim = "\t") %>%
  select(participant_id, germline_sample_platekey, tumour_sample_platekey)
dfMSS <- left_join(refitting, dfMSS) %>%
  filter(!is.na(germline_sample_platekey))

dfMSI <- read_delim("Timing_MSI_QCPASS.tsv", delim = "\t") %>%
  select(participant_id, germline_sample_platekey, tumour_sample_platekey)
dfMSI <- left_join(refitting, dfMSI) %>%
  filter(!is.na(germline_sample_platekey))


joined <- bind_rows(dfMSS, dfMSI) %>%
  select(participant_id, tumour_sample_platekey, germline_sample_platekey, karyotype)
write_delim(x = joined, "Timing_QCPASS_reruns.tsv", delim = "\t")
