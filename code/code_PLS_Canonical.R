#!/usr/bin/env Rscript
#This function executes the sPLS-Canonical and returns the redundancy
#The function takes as argument data as produced by code_simulation*.R and a microbiome normalization

args = commandArgs(trailingOnly=TRUE)

compute.redundancy.pls = function(X,Y,pls.object){
  
  
  X.scaled = apply(X, 2, function(col) (col-mean(col))/sd(col))
  Y.scaled = apply(Y, 2, function(col) (col-mean(col))/sd(col))
  
  min.col =  min(c(ncol(X),ncol(Y)))  
  #Proportion of variance by each variate
  variance.explained.X = sapply(1:min.col, function(y) sum(sapply(1:ncol(X), function(x) cor(X.scaled[,x], pls.object$variates$X[,y])^2))/ncol(X))
  variance.explained.Y = sapply(1:min.col, function(y) sum(sapply(1:ncol(Y), function(x) cor(Y.scaled[,x], pls.object$variates$Y[,y])^2))/ncol(Y))
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = sapply(1:min.col, function(x) cor(pls.object$variates$X[,x], pls.object$variates$Y[,x])^2)
  
  #Proportion of variance explained by the canonical variates
  cond.redundancy.X = variance.explained.X * proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y * proportion.variance.canonical.cor
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }


data = args[1]
microbiome.norm = args[2]

data = readRDS(data)
pls.results = lapply(1:1000, function(rep){
	print(rep)
	species = data[[rep]]$Simulated.Microbiotes
	metabolites = data[[rep]]$Simulated.Metabolites #should be log transformed
	
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

	min.col = min(c(ncol(normalized.species),ncol(metabolites)))

	pls = mixOmics::pls(normalized.species, metabolites, ncomp=min.col, mode="canonical")
	compute.redundancy.pls(normalized.species, metabolites, pls)
})

saveRDS(pls.results, paste0("results_PLS_canonical_",microbiome.norm,".RDS"))
