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
names(results_data) <- c("dataset_size","cdhit", "falconn")
results_data$cdhit <- results_data$cdhit / 60.0
results_data$falconn <- results_data$falconn / 60.0

p10 <- ggplot() + geom_line(data = results_data, aes(x = dataset_size, y = falconn, color="CSANN-Clust", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = falconn, color="CSANN-Clust", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = cdhit, color="CD-HIT", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = cdhit, color="CD-HIT", group="Algorithms"))
p10 <- p10 + ggtitle("") + labs(color = "Algorithms")
p10 <- p10 + xlab("Number of Database Sequences") + ylab("Clustering Time(M)")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_number_of_database_sequences_vs_clustering_time.png", sep=""), plot=p10)