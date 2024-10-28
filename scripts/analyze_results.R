library(zigzag)
library(dplyr)

setwd("~/Desktop/carrie/rhizome_trans")

# rhi 
nruns <- 6
prob_active_rhi <- vector("list", length = nruns)

for (i in 1:nruns) {
  print(i)
  path <- paste0("rhi", i, "/rhi", i, "_mcmc_output/rhi", i, "_probability_active.tab")
  prob_active_rhi[[i]] <- read.table(path,
                                     header = T, 
                                     row.names = 1)
}

# average across runs 
prob_active_avg_rhi <- (prob_active_rhi[[1]]+prob_active_rhi[[2]]+prob_active_rhi[[3]]+
                        prob_active_rhi[[4]]+prob_active_rhi[[5]]+prob_active_rhi[[6]])/6
colnames(prob_active_avg_rhi) <- "prob_active_rhi"
prob_active_avg_rhi$iso <- rownames(prob_active_avg_rhi)
rownames(prob_active_avg_rhi) <- NULL

# SAM
nruns <- 6
prob_active_SAM <- vector("list", length = nruns)

for (i in 1:nruns) {
  print(i)
  path <- paste0("SAM", i, "/SAM", i, "_mcmc_output/SAM", i, "_probability_active.tab")
  prob_active_SAM[[i]] <- read.table(path,
                                     header = T, 
                                     row.names = 1)
}

# average across runs 
prob_active_avg_SAM <- (prob_active_SAM[[1]]+prob_active_SAM[[2]]+prob_active_SAM[[3]]+
                          prob_active_SAM[[4]]+prob_active_SAM[[5]]+prob_active_SAM[[6]])/6
colnames(prob_active_avg_SAM) <- "prob_active_SAM"
prob_active_avg_SAM$iso <- rownames(prob_active_avg_SAM)
rownames(prob_active_avg_SAM) <- NULL

# tuber
nruns <- 6
prob_active_tub <- vector("list", length = nruns)

for (i in 1:nruns) {
  print(i)
  path <- paste0("tub", i, "/tub", i, "_mcmc_output/tub", i, "_probability_active.tab")
  prob_active_tub[[i]] <- read.table(path,
                                     header = T, 
                                     row.names = 1)
}

# average across runs 
prob_active_avg_tub <- (prob_active_tub[[1]]+prob_active_tub[[2]]+prob_active_tub[[3]]+
                        prob_active_tub[[4]]+prob_active_tub[[5]]+prob_active_tub[[6]])/6
colnames(prob_active_avg_tub) <- "prob_active_tub"
prob_active_avg_tub$iso <- rownames(prob_active_avg_tub)
rownames(prob_active_avg_tub) <- NULL


# root
nruns <- 6 #skip run 5 
prob_active_roo <- vector("list", length = nruns)

for (i in 1:nruns) {
  print(i)
  path <- paste0("roo", i, "/roo", i, "_mcmc_output/roo", i, "_probability_active.tab")
  prob_active_roo[[i]] <- read.table(path,
                                     header = T, 
                                     row.names = 1)
}

# average across runs 
prob_active_avg_roo <- (prob_active_roo[[1]]+prob_active_roo[[2]]+prob_active_roo[[3]]+
                        prob_active_roo[[4]]+prob_active_roo[[5]]+prob_active_roo[[6]])/6
colnames(prob_active_avg_roo) <- "prob_active_roo"
prob_active_avg_roo$iso <- rownames(prob_active_avg_roo)
rownames(prob_active_avg_roo) <- NULL

# combine datasets 
alpha1 <- 0.01
alpha5 <- 0.05
full_join(prob_active_avg_SAM, prob_active_avg_rhi) %>%
  full_join(prob_active_avg_roo) %>%
  full_join(prob_active_avg_tub) -> prob_active

prob_active$only_SAM <- prob_active$prob_active_SAM > 1-alpha1 & prob_active$prob_active_rhi < alpha1
prob_active$only_rhi <- prob_active$prob_active_rhi > 1-alpha1 & prob_active$prob_active_SAM < alpha1
prob_active$only_roo <- prob_active$prob_active_roo > 1-alpha5 & prob_active$prob_active_tub < alpha5
prob_active$only_tub <- prob_active$prob_active_tub > 1-alpha5 & prob_active$prob_active_roo < alpha5

# total number unique to SAM and rhi
sum(prob_active$only_SAM)
sum(prob_active$only_rhi)
sum(prob_active$only_roo)
sum(prob_active$only_tub)

# set of isoforms active in rhizomes and tubers (weird organ genes?)
USO_genes <- prob_active[prob_active$only_tub & prob_active$only_rhi, ]

# check annotations of these isoforms 
trans_ann <- read.csv("data/trinotate_annotation_report.csv")
write.csv(file = "USO_genes.csv", trans_ann[trans_ann$transcript_id %in% USO_genes$iso, ])
