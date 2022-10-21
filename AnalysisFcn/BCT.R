#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load libraries
library(BFpack)
library(R.matlab)

# Read data from MATLAB in version v6
tmpMat<-readMat(args[1])
df <- data.frame(tmpMat[[names(tmpMat)[1]]])

# Bayesian statistical tests
set.seed(123)
cor1 <- cor_test(df)

# Test specific hypothesis
if (length(args)>1) {	
	BF2 <- BF(cor1, hypothesis = args[2])
	# Print the results
	summary(BF2)
}	else {
	BF1 <- BF(cor1)
	summary(BF1)
} 