#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file file_prefix", call.=FALSE)
				}
				csv_file <- args[1]
				file_prefix <- args[2]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("dataset_size","PI70","PI80","PI90")

p10 <- ggplot() + geom_line(data = results_data, aes(x = dataset_size, y = PI70, color="> 70", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI70, color="> 70", group="threshold"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = PI80, color="> 80", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI80, color="> 80", group="threshold"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = PI90, color="> 90", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI90, color="> 90", group="threshold"))
p10 <- p10 + ggtitle("") + labs(color = "PI Threshold") +ylim(0,1)
p10 <- p10 + xlab("Number of Database Sequences") + ylab("Accuracy")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_number_of_database_sequences_vs_accuracy.png", sep=""), plot=p10)