#!/usr/bin/env Rscript
#This script executes the MMiRKAT method and produces the corresponding file for the microbiome normalization
#Different Distance kernels are considered
#Metabolome data can be log transformed See the corresponding line in the code

#The code takes 3 arguments 1) the microbiome normalization 2) the data as produced by simulate_data*.R scripts 3) scenario either null or alternative
library(MiRKAT)

apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }



args = commandArgs(trailingOnly=TRUE)

norm.Microbiome = args[1]
data = args[2]
scenario = args[3]


data = readRDS(data) #Should be a list as produced by simulate_data*.R scripts



results.MMiRKAT = lapply(1:length(data), function(rep) {
	print(rep)

	species = apply_nzv_filter(data[[rep]]$Simulated.Microbiotes)
	metabolites = apply_nzv_filter(data[[rep]]$Simulated.Metabolites) #Can be changed for considering the log metabolome

	if(norm.Microbiome == "clr"){
    		normalized.microbiome = compositions::clr(species)
  }
  	else if(norm.Microbiome =="ilr"){
    		normalized.microbiome = compositions::ilr(species)
  }
  	else if (norm.Microbiome == "alpha"){
    		best.alpha.microbiome = Compositional::alfa.tune(species + 1)[1]
    		normalized.microbiome = Compositional::alfa(species + 1, best.alpha.microbiome)$aff
  }

	dist.euclidean = as.matrix(dist(normalized.microbiome, method="euclidean"))
	dist.manhattan = as.matrix(dist(normalized.microbiome, method="manhattan"))
	dist.canberra = as.matrix(dist(normalized.microbiome, method="canberra"))
	
	dist.euclidean[is.na(dist.euclidean)] = 0
	dist.manhattan[is.na(dist.manhattan)] = 0
	dist.canberra[is.na(dist.canberra)] = 0

	KS = list("euclidean"=D2K(dist.euclidean), "manhattan"=D2K(dist.manhattan), "canberra"=D2K(dist.canberra))
	MMiRKAT(Y=metabolites,Ks=KS)$p_values
})

saveRDS(results.MMiRKAT, paste0("results_MMiRKAT_",norm.Microbiome,"_",scenario, ".RDS"))
	






