library(Biobase)
require(foreign)
require(MASS)
library(data.table)

print('start')

#oarse snakemake input/params
signal_all <- snakemake@input[["signal"]]
regVars <- snakemake@input[["covariates"]]
include <- snakemake@input[["include"]]
outSignal <- snakemake@output[["outSignal"]]
marks <- snakemake@params[["marks"]]
chrm <- snakemake@params[["chrm"]]
#variables to regress out 
regVars <- read.table(regVars, header=T)

#center and scale numerical variables
regVars$numReads <- scale(regVars$numReads, center=T, scale=T)
regVars$RSC <- scale(regVars$RSC, center=T, scale=T)

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
