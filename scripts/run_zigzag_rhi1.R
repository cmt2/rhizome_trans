#zigzag analysis of diff expr. for rhizome vs. sam 
library(zigzag)

setwd("~/Desktop/carrie/rhizome_trans")

# set up run info: 
r <- 1 # r <- 2, r <- 3, r <- 4
tissue <- "rhi"# tissue <- "SAM", tissue <- "roo", tissue <- "tub
print(paste0("Running zigzag on ", tissue, ", run ", r))

txi.rsem <- readRDS("data/processed_data_txi.RDS")

# get in zigzag format: 
quant <- txi.rsem$abundance[, grep(tissue, x = colnames(txi.rsem$abundance))]
lengths <- data.frame(length = rowMeans(txi.rsem$length))

print("setting up object")
myobj <- zigzag$new(quant, 
                    lengths,
                    threshold_a = c(0,3.5), 
                    active_variances_prior_max = 3.5,
                    output_directory = paste0(tissue, r))

print("running burnin")
myobj$burnin(sample_frequency = 50, 
             ngen=20000, 
             write_to_files = TRUE)

print("running MCMC")
myobj$mcmc(sample_frequency = 50, 
           ngen = 200000, 
           run_posterior_predictive = TRUE, 
           mcmcprefix = paste0(tissue, r))