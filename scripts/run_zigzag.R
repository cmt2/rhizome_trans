#zigzag analysis of diff expr. for rhizome vs. sam 
library(tximport)
library(tidyverse)
library(ggthemes)
library(zigzag)

setwd("~/Documents/rhizome_trans")

conditions <- read.table("data/samples_described.txt",sep="\t", header = TRUE,
                         row.names = 1)
files <- file.path(paste0("data/isoform_counts/RSEM.isoforms.results.", 
                   rownames(conditions)))
names(files) <- rownames(conditions)

txi.rsem <- tximport(files, 
                     txIn = TRUE, 
                     txOut = TRUE,  
                     type = "rsem")
saveRDS(txi.rsem, file = "data/processed_data_txi.RDS")
readRDS("data/processed_data_txi.RDS")
# get in zigzag format: 

sam_quant <- txi.rsem$abundance[,1:3]
rhi_quant <- txi.rsem$abundance[,4:6]
lengths <- data.frame(length = rowMeans(txi.rsem$length))

sam_zigzag <- zigzag$new(sam_quant, lengths)
#sam_zigzag_run2 <- zigzag$new(sam_quant, lengths)
#rhi_zigzag_run1 <- zigzag$new(rhi_quant, lengths)
#rhi_zigzag_run2 <- zigzag$new(rhi_quant, lengths)

sam_zigzag$burnin(sample_frequency = 50, ngen=20000, write_to_files = TRUE)
sam_zigzag$mcmc(sample_frequency = 50, ngen = 50000, run_posterior_predictive = TRUE, mcmcprefix = "SAM_run1")


