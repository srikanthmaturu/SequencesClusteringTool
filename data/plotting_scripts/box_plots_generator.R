#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file plot_title.n", call.=FALSE)
				}
				csv_file <- args[1]
				plot_title <- args[2]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("edit_distance", "percent_similarity", "pure_cosine_distance", "cosine_distance", "pure_cosine_angle", "cosine_angle", "pure_euclidean_distance", "euclidean_distance")

#cosine_data <- unique(results_data[,c("edit_distance","cosine_distance")])

#cosine_data <- cosine_data[order(cosine_data$edit_distance),]

#png(file = paste(plot_title, "_Cosine_Distance.png", sep=""))
#boxplot(cosine_distance ~ edit_distance, data = cosine_data, ylab = "Cosine distance", xlab = "Edit distance", main=plot_title)
#dev.off()

#rm(cosine_data)

#euc_data <- unique(results_data[,c("edit_distance","euclidean_distance")])

#euc_data <- euc_data[order(euc_data$edit_distance),]

#png(file = paste(plot_title, "_Euclidean_Distance.png", sep=""))
#boxplot(euclidean_distance ~ edit_distance, data = euc_data, ylab = "Euclidean distance", xlab = "Edit distance", main=plot_title)
#dev.off()

#rm(euc_data)
#abline(h=mean(x))

percent_sim_data <- unique(results_data[,c("percent_similarity","cosine_distance")])

percent_sim_data <- percent_sim_data[order(percent_sim_data$percent_similarity),]

png(file = paste(plot_title, "_Cosine_Distance_vs_Percent_Similarity.png", sep=""))
boxplot(cosine_distance ~ percent_similarity, data = percent_sim_data, ylab = "Cosine Similarity", xlab = "Percent Identity", main=plot_title)
means <- tapply(percent_sim_data$cosine_distance,percent_sim_data$percent_similarity,mean)
points(means,col="red",pch=18)
dev.off()

rm(percent_sim_data)

euc_percent_sim_data <- unique(results_data[,c("percent_similarity","euclidean_distance")])

euc_percent_sim_data <- euc_percent_sim_data[order(euc_percent_sim_data$percent_similarity),]

png(file = paste(plot_title, "_Euclidean_Distance_vs_Percent_Similarity.png", sep=""))
boxplot(euclidean_distance ~ percent_similarity, data = euc_percent_sim_data, ylab = "Squared Euclidean Distance", xlab = "Percent Identity", main=plot_title)
means <- tapply(euc_percent_sim_data$euclidean_distance,euc_percent_sim_data$percent_similarity,mean)
points(means,col="red",pch=18)
dev.off()


q(status=0)