#!/usr/bin/env Rscript
#This function executes the procrustes analysis 
#This function takes as argument data as produced by code_simulation*.R and a microbiome normalization
args = commandArgs(trailingOnly=TRUE)


library(vegan)


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
  nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
  if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
  return(x)
}


data = readRDS(args[1])
norm.Microbiome = args[2]

procrustes.results = sapply(1:length(data), function(x){
	
	print(x)

	species = data[[x]]$Simulated.Microbiotes
	metabolites = data[[x]]$Simulated.Metabolites #should be log transformed
	
	if(norm.Microbiome == "clr"){
    		normalized.microbiome = compositions::clr(Microbiome)
  	} 
  	else if(norm.Microbiome =="ilr"){
    		normalized.microbiome = compositions::ilr(Microbiome)
  	} 
  	else if (norm.Microbiome == "alpha"){
    		best.alpha.microbiome = Compositional::alfa.tune(Microbiome + 1)[1]
    		normalized.microbiome = Compositional::alfa(Microbiome + 1, best.alpha.microbiome)$aff
  	}	
		
	normalized.microbiome = apply_nzv_filter(normalized.microbiome)
	metabolites = apply_nzv_filter(metabolites)

	species.rda = rda(normalized.microbiome)
	metabolites.rda = rda(metabolites)
	
	protest(X = species.rda, Y = metabolites.rda)$signif
	})

saveRDS(procrustes.results, paste0("Procrustes_", norm.Microbiome ,".RDS"))
