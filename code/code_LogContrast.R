#!/usr/bin/env Rscript
#This script executes the log contrast regression and produces the corresponding p-values for the ANOVA test
#The script takes data as produced by code_simulation*.R as argument

args = commandArgs(trailingOnly=TRUE)


set.seed(1234)
library(Compositional)

compute_F_from_LCLR = function(y, model){
  mean.y = mean(y, na.rm=T)
  ss.total = sum((y - mean.y)^2)
  ss.residuals = sum(model$residuals^2)
  ss.model = ss.total - ss.residuals
  ms.model = ss.model / (length(model$be)-1)
  ms.residuals = ss.residuals/(length(y)-length(model$be))
  F_stat = ms.model/ms.residuals

  return(pf(F_stat,(length(model$be)-1), (length(y)-length(model$be)),lower.tail = F))


}

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }

data = readRDS(args[1])

log.contrast.results = lapply(1:length(data), function(x) {
	print(x)
	species = as.matrix(data[[x]]$Simulated.Microbiotes+1)
	metabolites = as.matrix(data[[x]]$Simulated.Metabolites) #should be log-transformed in some scenarios
	
	colnames(species) = paste0("Micro_", 1:ncol(species))
	colnames(metabolites) = paste0("Metabo_", 1:ncol(metabolites))
	
	species = apply_nzv_filter(species)
	metabolites = apply_nzv_filter(metabolites)

	list.results = list()

	for(col in colnames(metabolites)){
        	
        	model = lc.reg(y=metabolites[,col], x=species)
        	pvalue = compute_F_from_LCLR(metabolites[,col], model)

        	list.results[[col]] = list("pvalue"=pvalue)
	}
	
	list.results



})


saveRDS(log.contrast.results,"results_logcontrast.RDS")
