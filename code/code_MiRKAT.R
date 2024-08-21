#!/usr/bin/env Rscript

#This script executes the MiRKAT method and produces the corresponding file for the microbiome normalization
#Different Distance kernels are considered
#Metabolome data can be log transformed See the corresponding line in the code

#The code takes 3 arguments 1) the microbiome normalization 2) the data as produced by simulate_data*.R scripts 3) scenario either null or alternative

args = commandArgs(trailingOnly=TRUE)
set.seed(1234)

library(MiRKAT)


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }

data = readRDS(args[1])
norm.Microbiome = args[2]
scenario = args[3]

MiRKAT.results = lapply(1:length(data), function(x) {
	print(x)

	species = data[[rep]]$Simulated.Microbiotes

	if(norm.Microbiome == "clr"){
    		species = compositions::clr(species)
  }
  	else if(norm.Microbiome =="ilr"){
    		species = compositions::ilr(species)
  }
  	else if (norm.Microbiome == "alpha"){
    		best.alpha.microbiome = Compositional::alfa.tune(species + 1)[1]
    		species = Compositional::alfa(species + 1, best.alpha.microbiome)$aff
  }

	metabolites = as.matrix(data[[x]]$Simulated.Metabolites)
	
	colnames(species) = paste0("Micro_", 1:ncol(species))
	colnames(metabolites) = paste0("Metabo_", 1:ncol(metabolites))
	
	species = apply_nzv_filter(species)
	metabolites = apply_nzv_filter(metabolites)

	list.results = list()
	
	dist.euclidean = as.matrix(dist(species, method="euclidean"))
	dist.manhattan = as.matrix(dist(species, method="manhattan"))
	dist.canberra = as.matrix(dist(species, method="canberra"))
	
	dist.euclidean[is.nan(dist.euclidean)] = 0
	dist.manhattan[is.nan(dist.manhattan)] = 0
	dist.canberra[is.nan(dist.canberra)] = 0
	
	K.list = list("E"=D2K(dist.euclidean),"M"=D2K(dist.manhattan),"C"=D2K(dist.canberra))


	for(col in colnames(metabolites)){
        	

        	list.results[[col]] = MiRKAT(y=metabolites[,col], Ks=K.list, out_type="C", omnibus="cauchy")
	}
	
	list.results



})


saveRDS(MiRKAT.results,"results_MiRKAT_", scenario,"_" ,norm.Microbiome,".RDS")
