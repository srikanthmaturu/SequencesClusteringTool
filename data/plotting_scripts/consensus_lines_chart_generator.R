#! /usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
if(length(args)<2) {
                stop("Usage ./executable.R csv_file file_prefix", call.=FALSE)
				}
				csv_file <- args[1]
				true_positives <- as.integer(args[2])
				file_prefix <- args[3]

results_data <- read.csv(csv_file, header=FALSE)
names(results_data) <- c("percent_identity", "Trinity",	"IDBA",	"SOAPdenovo",	"SPAdes",	"consensus")
results_data <- results_data[order(-as.numeric(results_data$percent_identity)),]
results_data$Trinity <- cumsum(results_data$Trinity)
results_data$IDBA <- cumsum(results_data$IDBA)
results_data$SOAPdenovo <- cumsum(results_data$SOAPdenovo)
results_data$SPAdes <- cumsum(results_data$SPAdes)
results_data$consensus <- cumsum(results_data$consensus)

p10 <- ggplot() + geom_line(data = results_data, aes(x = percent_identity, y = Trinity, color="Trinity", group="Assemblers")) + geom_point(data = results_data, aes(x = percent_identity, y = Trinity, color="Trinity", group="Assemblers"))
p10 <- p10 + geom_line(data = results_data, aes(x = percent_identity, y = IDBA, color="IDBA", group="Assemblers")) + geom_point(data = results_data, aes(x = percent_identity, y = IDBA, color="IDBA", group="Assemblers"))
p10 <- p10 + geom_line(data = results_data, aes(x = percent_identity, y = SOAPdenovo, color="SOAPdenovo", group="Assemblers")) + geom_point(data = results_data, aes(x = percent_identity, y = SOAPdenovo, color="SOAPdenovo", group="Assemblers"))
p10 <- p10 + geom_line(data = results_data, aes(x = percent_identity, y = SPAdes, color="SPAdes", group="Assemblers")) + geom_point(data = results_data, aes(x = percent_identity, y = SPAdes, color="SPAdes", group="Assemblers"))
p10 <- p10 + geom_line(data = results_data, aes(x = percent_identity, y = consensus, color="All Combined", group="Assemblers")) + geom_point(data = results_data, aes(x = percent_identity, y = consensus, color="All Combined", group="Assemblers"))
p10 <- p10 + geom_hline(yintercept = true_positives, color="black") + geom_text(aes(80,true_positives,label=paste("Total True Positives = ", true_positives, sep=""),vjust=-1))
p10 <- p10 + ggtitle("") + xlim(100,70) + labs(color = "Assemblers")
p10 <- p10 + scale_y_continuous(sec.axis = sec_axis(~./(as.integer(true_positives/100)), name = "Cumulative Percentage of True Positives Found [%]"))
p10 <- p10 + xlab("Percent Identity") + ylab("Cumulative Number of True Positives Found")
p10 <- p10 + theme_classic()
ggsave(paste(file_prefix, "_consensus_percent_identity_vs_total_true_positives_found.png", sep=""), plot=p10)