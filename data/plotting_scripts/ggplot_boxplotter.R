#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<1) {
                stop("Usage ./executable.R csv_file", call.=FALSE)
				}
				csv_file <- args[1]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("edit_distance", "percent_similarity", "pure_cosine_distance", "cosine_distance", "pure_cosine_angle", "cosine_angle", "pure_euclidean_distance", "euclidean_distance")

percent_sim_data <- unique(results_data[,c("percent_similarity","cosine_distance")])
percent_sim_data <- percent_sim_data[order(percent_sim_data$percent_similarity),]

p10 <- ggplot(aes( x = factor(percent_similarity), y = cosine_distance, group = percent_similarity), data = percent_sim_data) + geom_boxplot()
p10 <- p10 + scale_x_discrete(name = "Percent Identity", breaks = seq(0,100,5)) + scale_y_continuous(name = "Cosine Similarity",
                              breaks = seq(-0.2, 1.0, 0.1),
                              limits=c(-0.2, 1.0))
p10 <- p10 + ggtitle("")
p10 <- p10 + stat_summary(fun.y=mean, geom="line", color="red", aes(group=1))
p10 <- p10 + theme_classic()
ggsave("percent_identity_vs_cosine_similarity.png", plot=p10)

rm(percent_sim_data)

euc_percent_sim_data <- unique(results_data[,c("percent_similarity","euclidean_distance")])
euc_percent_sim_data <- euc_percent_sim_data[order(euc_percent_sim_data$percent_similarity),]

p10 <- ggplot(aes( x = factor(percent_similarity), y = euclidean_distance, group = percent_similarity), data = euc_percent_sim_data) + geom_boxplot()
p10 <- p10 + scale_x_discrete(name = "Percent Identity", breaks = seq(0,100,5)) + scale_y_continuous(name = "Squared Euclidean Distance",
                              breaks = seq(0, 2, 0.1),
                              limits=c(0, 2))
p10 <- p10 + ggtitle("")
p10 <- p10 + stat_summary(fun.y=mean, geom="line", color="red", aes(group=1))
p10 <- p10 + theme_classic()
ggsave("percent_identity_vs_squared_euclidean_distance.png", plot=p10)
