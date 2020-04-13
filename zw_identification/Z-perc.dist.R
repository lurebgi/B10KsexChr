args <- commandArgs(trailingOnly=TRUE)

library(ggplot2)

data <- read.table(args[1])

pdf(args[2], height=6,width=8)
ggplot(data) + geom_histogram(aes(x=V3),binwidth=1)

dev.off()
