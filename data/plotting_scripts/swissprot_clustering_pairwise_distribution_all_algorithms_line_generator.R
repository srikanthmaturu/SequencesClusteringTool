#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
print(args)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file file_prefix", call.=FALSE)
				}
				csv_file <- args[1]
				file_prefix <- args[2]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("id","cluster_size_csann","min_pi_csann","average_pi_csann", "cluster_size_cdhit","min_pi_cdhit","average_pi_cdhit")
results_data$avg_pairwise_distance_csann <- 1 - results_data$average_pi_csann / 100.0
results_data$avg_pairwise_distance_cdhit <- 1 - results_data$average_pi_cdhit / 100.0


p10 <- ggplot() + geom_line(data = results_data, aes(x = avg_pairwise_distance_csann, y = cluster_size_csann, color="CSANN-Clust", group="Algorithms")) + geom_point(data = results_data, aes(x = avg_pairwise_distance_csann, y = cluster_size_csann, color="CSANN-Clust", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = avg_pairwise_distance_cdhit, y = cluster_size_cdhit, color="CD-HIT", group="Algorithms")) + geom_point(data = results_data, aes(x = avg_pairwise_distance_cdhit, y = cluster_size_cdhit, color="CD-HIT", group="Algorithms"))
p10 <- p10 + ggtitle("") + labs(color = "Algorithms")
p10 <- p10 + xlab("Average pairwise distance") + ylab("Relative Frequency")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_relative_frequency_vs_average_pairwise_distance.png", sep=""), plot=p10)