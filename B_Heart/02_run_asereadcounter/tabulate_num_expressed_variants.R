#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1], sep = "\t")

data_subset = subset(data, data$totalCount>10)

out = data.frame(counts = c(nrow(data), nrow(data_subset)), labels = c("before_filtering", "after_filtering"))
# write.table(out, args[2], quote = F, sep = "\t", row.names = F)
print(out)