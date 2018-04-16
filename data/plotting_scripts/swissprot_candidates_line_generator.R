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
names(results_data) <- c("dataset_size","PI70","PI80","PI90")
results_data$PI70 = results_data$PI70/query_size
results_data$PI80 = results_data$PI80/query_size
results_data$PI90 = results_data$PI90/query_size

p10 <- ggplot() + geom_line(data = results_data, aes(x = dataset_size, y = PI70, color="> 70", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI70, color="> 70", group="threshold"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = PI80, color="> 80", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI80, color="> 80", group="threshold"))
p10 <- p10 + geom_line(data = results_data, aes(x = dataset_size, y = PI90, color="> 90", group="threshold")) + geom_point(data = results_data, aes(x = dataset_size, y = PI90, color="> 90", group="threshold"))
p10 <- p10 + ggtitle("") + labs(color = "PI Threshold")
p10 <- p10 + xlab("Number of Database Sequences") + ylab("Average No. of Candidates Sequences")
#p10 <- p10 + scale_y_continuous(sec.axis = sec_axis(~./(as.integer(results_data$dataset_size))), name = "No. of Candidates/Database Size")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_number_of_database_sequences_vs_number_of_candidates.png", sep=""), plot=p10)
