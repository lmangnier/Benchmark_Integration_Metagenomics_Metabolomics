#!/usr/bin/env Rscript
#This function executes the RDA and returns the explained variance 
#The function takes as argument data as provided by code_simulation*.R and a microbiome normalization
args = commandArgs(trailingOnly=TRUE)
library(vegan)

data = readRDS(args[1])
microbiome.norm = args[2]


compute.explained.variance.rda = function(rda.object){
  
  return(rda.object$CCA$tot.chi/ rda.object$tot.chi)
}

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }

results.RDA = lapply(1:length(data), function(rep) {
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

	rda = vegan::rda(normalized.species, metabolites) #species and metabolites can be switched if needed
  	compute.explained.variance.rda(rda)
	
})

saveRDS(results.RDA, paste0("RDA_",microbiome.norm,".RDS"))
