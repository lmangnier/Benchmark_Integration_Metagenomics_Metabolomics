#!/usr/bin/env Rscript
#This function returns the results for the CLR-MLASSO
#The function takes as an argument data as provided by the code_simulation*.R

library(glmnet)


args = commandArgs(trailingOnly=TRUE)

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
  nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
  if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
  return(x)}

data = args[1]

data = readRDS(data)

clr.MLASSO.results = lapply(1:length(data), function(rep){
        print(rep)
        species = compositions::clr(data[[rep]]$Simulated.Microbiotes)
        metabolites = data[[rep]]$Simulated.Metabolites #log transformation should be considered

        colnames(species) = paste0("Micro", 1:ncol(species))
        colnames(metabolites) = paste0("Metabo", 1:ncol(metabolites))

        species = apply_nzv_filter(species)
        metabolites = apply_nzv_filter(metabolites)
	
	cvfit = cv.glmnet(species, metabolites,family="mgaussian", alpha=1, type.measure="mse")
        coefs.clr.MLASSO=coef(cvfit, s="lambda.min")

	return(coefs.clr.MLASSO)
})

saveRDS(clr.MLASSO.results, "results_clr_MLASSO.RDS")

