#!/usr/bin/Rscript
#Perform statistical tests for varstats.py
#wdecoster

args <- commandArgs(trailingOnly = TRUE)
data <- matrix(sapply(strsplit(args, " "), as.numeric), nrow = 2, byrow=TRUE)

cat(fisher.test(data)['p.value']$p.value)
#data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow=TRUE)
cat(" ")
suppressMessages(library(HardyWeinberg))
suppressWarnings(cat(HWExact(colSums(data), verbose=F)$pval))