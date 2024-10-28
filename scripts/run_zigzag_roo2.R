#zigzag analysis of diff expr. for rhizome vs. sam 
library(zigzag)

setwd("~/Desktop/carrie/rhizome_trans")

# set up run info: 
r <- 2 # r <- 1, r <- 3, r <- 4
tissue <- "roo"# tissue <- "rhi", tissue <- "SAM", tissue <- "tub
print(paste0("Running zigzag on ", tissue, ", run ", r))

txi.rsem <- readRDS("data/processed_data_txi.RDS")

# get in zigzag format: 
quant <- txi.rsem$abundance[, grep(tissue, x = colnames(txi.rsem$abundance))]
lengths <- data.frame(length = rowMeans(txi.rsem$length))

print("setting up object")
myobj <- zigzag$new(quant, 
                    lengths,
                    active_variances_prior_max = 2.5,
                    threshold_a = c(1, 3),
                    output_directory = paste0(tissue, r))

print("running burnin")
myobj$burnin(sample_frequency = 50, 
             ngen=5000, # reduced for roo samples 
             write_to_files = TRUE)

print("running MCMC")
myobj$mcmc(sample_frequency = 50, 
           ngen = 200000, 
           run_posterior_predictive = TRUE, 
           mcmcprefix = paste0(tissue, r))