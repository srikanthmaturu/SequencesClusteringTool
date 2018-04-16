#! /usr/bin/env Rscript
library(ggplot2)

##install.packages("ggpubr")
library(ggpubr)


args <- commandArgs(TRUE)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file plot_title.n", call.=FALSE)
				}
				csv_file <- args[1]
				plot_title <- args[2]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("edit_distance", "percent_similarity", "pure_cosine_distance", "cosine_distance", "pure_cosine_angle", "cosine_angle", "pure_euclidean_distance", "euclidean_distance")


percent_sim_data <- unique(results_data[,c("percent_similarity","cosine_distance")])

percent_sim_data <- percent_sim_data[order(percent_sim_data$percent_similarity),]

percent_sim_data <- percent_sim_data[percent_sim_data$percent_similarity>=50,]

#cor.test(percent_sim_data$percent_similarity,percent_sim_data$cosine_distance, method=c("pearson", "kendall", "spearman"))
print("cosine distance vs percent similarity")
cor.test(percent_sim_data$cosine_distance, percent_sim_data$percent_similarity, method=c("pearson"))
cor.test(percent_sim_data$cosine_distance, percent_sim_data$percent_similarity, method=c("kendall"))
cor.test(percent_sim_data$cosine_distance, percent_sim_data$percent_similarity, method=c("spearman"))



rm(percent_sim_data)

euc_percent_sim_data <- unique(results_data[,c("percent_similarity","euclidean_distance")])

euc_percent_sim_data <- euc_percent_sim_data[order(euc_percent_sim_data$percent_similarity),]

euc_percent_sim_data <- euc_percent_sim_data[euc_percent_sim_data$percent_similarity>=50,]

print("Euclidean vs percent similarity")
cor.test(euc_percent_sim_data$euclidean_distance,euc_percent_sim_data$percent_similarity, method=c("pearson"))
cor.test(euc_percent_sim_data$euclidean_distance,euc_percent_sim_data$percent_similarity, method=c("kendall"))
cor.test(euc_percent_sim_data$euclidean_distance,euc_percent_sim_data$percent_similarity, method=c("spearman"))


#cor.test(euc_percent_sim_data$percent_similarity, euc_percent_sim_data$euclidean_distance, method=c("pearson","kendall","spearman"))

q(status=0)