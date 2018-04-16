#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<3) {
                stop("Usage ./executable.R csv_file queries_size file_prefix", call.=FALSE)
				}
				csv_file <- args[1]
				query_size <- as.integer(args[2])
				file_prefix <- args[3]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("dataset_size", "blast", "edlib", "falconn_PI70", "falconn_PI80", "falconn_PI90")
results_data$falconn_PI70 = results_data$falconn_PI70/query_size
results_data$falconn_PI80 = results_data$falconn_PI80/query_size
results_data$falconn_PI90 = results_data$falconn_PI90/query_size
results_data$edlib = results_data$edlib/query_size
results_data$blast = results_data$blast/query_size

p10 <- ggplot() + geom_line(data = results_data, aes(x = dataset_size, y = falconn_PI70, color="CSANN(PI > 70)", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = falconn_PI70, color="CSANN(PI > 70)", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = falconn_PI80, color="CSANN(PI > 80)", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = falconn_PI80, color="CSANN(PI > 80)", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = falconn_PI90, color="CSANN(PI > 90)", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = falconn_PI90, color="CSANN(PI > 90)", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = edlib, color="Brute-force(NW)", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = edlib, color="Brute-force(NW)", group="Algorithms"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = blast, color="BLASTP", group="Algorithms")) + geom_point(data = results_data, aes(x = dataset_size, y = blast, color="BLASTP", group="Algorithms"))
p10 <- p10 + ggtitle("") + labs(color = "Algorithms")
p10 <- p10 + xlab("Number of Database Sequences") + ylab("Query Time (sec)")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_number_of_database_sequences_vs_query_time.png", sep=""), plot=p10)