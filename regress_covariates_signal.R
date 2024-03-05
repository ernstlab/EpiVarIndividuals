library(Biobase)
require(foreign)
require(MASS)
library(data.table)


args = commandArgs(trailingOnly=TRUE)
signal_all = args[1] #Signal output from ChromHMM, must have mark in sample names
regVars = args[2] #Table of covariates 
outSignal = args[3] #output file for corrected signal
marks = args[4] #Comma separated string of marks
chrm = args[5] #Chromosome

regVars <- read.table(regVars, header=T)
marks <- strsplit(marks, ",")[[1]]

#center and scale numerical variables
ind_scale <- sapply(regVars, is.numeric)
regVars[ind_scale] <- lapply(regVars[ind_scale], scale)

#uncorrected signal
signal_all <- read.table(signal_all, header=T, skip=1)
print(dim(signal_all))
signal_out <- list()
samples_all <- c()
for(mark in marks){
	print(mark)

	#subset signal to include only one mark
	keep <- c()
	for(c in colnames(signal_all)){
		if (grepl( mark, c, fixed = TRUE)){
			keep <- c(keep, c)
		}
	}
	signal <- signal_all[,keep]

	
	#get overlapping samples
	samples <- intersect(rownames(regVars), colnames(signal))	
	samples_all <- c(samples_all, samples)
	
	signal <- signal[,samples]
	regVars_sub <- regVars[samples,]
	
	#quantile normalize
	signal <- data.matrix(signal)
	quantile_normalisation <- function(df){
	  df_rank <- apply(df,2,rank,ties.method="min")
	  df_sorted <- data.frame(apply(df, 2, sort))
	  df_mean <- apply(df_sorted, 1, mean)
	   
	  index_to_mean <- function(my_index, my_mean){
	    return(my_mean[my_index])
	  }
	   
	  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
	  rownames(df_final) <- rownames(df)
	
	  return(df_final)
	}
	
	signal <- quantile_normalisation(signal)
	
	#round to make count data
	signal <- round(signal)
	#peudo count
	signal <- signal + 1
	
	#Regress out covariates
	# Turn warnings into errors so they can be trapped
	options(warn = 2)
	control0 <- glm.control(epsilon = 1e10, maxit = 100) #increase the default number of iterations from 25 to 100
	adjust <- function(bin){
	        fit <- glm(bin ~ . , data = regVars_sub, family=quasipoisson)
	        return((bin*exp(coefficients(fit)[1]))/predict(fit, type="response"))
	}
	
	signalCorrected <- round(apply(signal, 1, adjust))
	signalCorrected <- data.frame(signalCorrected)
	rownames(signalCorrected) = colnames(signal)
	signal_out[[mark]] = signalCorrected
	
}

signal_out_all <- rbindlist(signal_out)
signal_out_all <- t(signal_out_all)
print(dim(signal_out_all))

#write corrected signal
header <- paste0("cell1\t", "chr", chrm)
cols <- paste(samples_all, sep="\t", collapse="\t")
#cat(header, file = outSignal)
writeLines(c(header, cols), outSignal, sep="\n")
#writeLines(samples_all, outSignal, sep="\t")
write.table(signal_out_all, file=outSignal, row.names=F, col.names=F, quote=F, sep="\t", append = TRUE)
