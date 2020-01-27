Rlib <- "~/re_gecip/cancer_colorectal/analysisResults/6.evolutionaryAnalyses/A.ChromosomeTiming/BATCH/code/Rpackages_3.5"
library(dplyr, lib.loc = Rlib)
library(ggplot2, lib.loc = Rlib)
library(tidyr, lib.loc = Rlib)
library(ggfortify, lib.loc = Rlib)
library(bbmle, lib.loc = Rlib)
library(ggpubr, lib.loc = Rlib)
library(sads, lib.loc = Rlib)
library(ggthemes, lib.loc = Rlib)
library(pio, lib.loc = Rlib)
library(easypar, lib.loc = Rlib)
library(mobster, lib.loc = Rlib)

basedir <- "~/re_gecip/shared_allGeCIPs/R_dm/data/COLORECTAL_timing/"
load("Timing_Status.RData")

df <- tibble()

for (i in 1:length(status)){
  if (status[[i]]$Timing_plot == FALSE){
    message("No timing results")
    next
  }
  if (round(i / 25) == (i / 25)){
    message(paste0(i, "/", length(status), " samples finished"))
  }
  id <- status[[i]]$participant_id
  fits <- list.files(paste0(basedir, id), pattern = "fit_and_data")
  for (f in fits){
    load(paste0(basedir, id, "/", f))
    karyotype <- str_split(f, "_")[[1]][5]
    print(paste0("ID = ", id, ", Karyotype = ", karyotype))
    newfit <- mobster:::rename_Beta_clusters(plot$fit$best)
    dftemp <- mobster:::to_string(newfit) %>%
      mutate(id = id,
             karyotype = karyotype)
    dftemp <- bind_cols(dftemp, plot$fit$best$scores)
    df <- bind_rows(df, dftemp)
  }
}

readr::write_csv(df, "timing_QC_table.csv")
