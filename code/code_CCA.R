#!/usr/bin/env Rscript
#This code executes the CCA producing the redundancy index for a given microbiome normalization
#The script takes as argument the data as produced by code_simulation*.R and a microbiome normalization

args = commandArgs(trailingOnly=TRUE)

library(CCA)

compute.redundancy.cca = function(X,Y,cca.object){
  
  #Proportion of variance by each variate
  variance.explained.X = colSums(cca.object$scores$corr.X.xscores^2)/ncol(X)
  variance.explained.Y = colSums(cca.object$scores$corr.Y.yscores^2)/ncol(Y)
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = cca.object$cor^2
  
  #Proportion of X explained by the canonical variates 
  cond.redundancy.X = variance.explained.X*proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y*proportion.variance.canonical.cor
  
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){   
	nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)   
	if(length(nzv$Position) > 0) x <- x[, -nzv$Position]   
	return(x) }


data = args[1]
microbiome.norm = args[2]

data = readRDS(data)
cca.results = lapply(1:length(data), function(rep){
	print(rep)
	species = data[[rep]]$Simulated.Microbiotes
	metabolites = data[[rep]]$Simulated.Metabolites #should be log-transformed in some scenarios
	
	colnames(species) = paste0("Micro", 1:ncol(species))
	colnames(metabolites) = paste0("Metabo", 1:ncol(metabolites))	
	
	species = apply_nzv_filter(species)
	metabolites = apply_nzv_filter(metabolites)
	
	
	if(microbiome.norm=="clr") normalized.species = compositions::clr(species)
	if(microbiome.norm=="ilr") normalized.species = compositions::ilr(species)
	if(microbiome.norm=="alpha"){
		 
		 best.alpha.microbiome = Compositional::alfa.tune(species + 1)[1]
                 normalized.species = Compositional::alfa(species + 1, best.alpha.microbiome)$aff
	}


	cc = tryCatch(cc(normalized.species, metabolites),  error=function(err) NA) 
	list("Redundancy" = tryCatch(compute.redundancy.cca(normalized.species, metabolites, cc), error=function(err) NA), "n.index.kept" = list("Microbiome"=ncol(species),"Metabolome"=ncol(metabolites)))
})

saveRDS(cca.results, paste0("results_CCA_",microbiome.norm,".RDS"))
