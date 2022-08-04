#!/usr/bin/env R

# arg1: R memory image of trained SOM
# arg2: catalogue

library(kohonen)
library(data.table)
library(helpRfuncs)


inputs<-commandArgs(T)

somfile <- inputs[1]
cat(paste("=> loading SOM file:", somfile, "\n"))
load(somfile)

datafile <- inputs[2]
cat(paste("=> loading catalogue data:", datafile, "\n"))
phot <- helpRfuncs::read.file(datafile)
som_wphot <- kohparse(train.som, phot, n.cores=64)
# write out the cell index
outend <- helpRfuncs::vecsplit(datafile, ".", -1)
outfile <- sub(paste0(".", outend), "_cellIDX.txt", datafile)
cat(paste("=> writing cell indices to ASCII file:", outfile, "\n"))
cat(sapply(som_wphot$unit.classif, toString), file=outfile, sep="\n")
