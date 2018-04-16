#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file file_prefix", call.=FALSE)
				}
				csv_file <- args[1]
				file_prefix <- args[2]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("percent_identity", "edlib_results", "falconn_results")
results_data$accuracy = results_data$falconn_results/results_data$edlib_results
total_falconn_items = sum(results_data$falconn_results)
total_edlib_items = sum(results_data$edlib_results)
overall_accuracy = total_falconn_items / total_edlib_items

p10 <- ggplot(data = results_data, aes(x = percent_identity, y = accuracy)) + geom_line() + geom_point()
p10 <- p10 + geom_hline(yintercept = overall_accuracy, color="red") + geom_text(aes(90,overall_accuracy,label=paste("Overall accuracy = ", round(overall_accuracy,2), sep=""),vjust=-1))
p10 <- p10 + ggtitle("")
p10 <- p10 + xlab("Percent Identity") + ylab("Accuracy")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_accuracy_for_percent_identity.png", sep=""), plot=p10)