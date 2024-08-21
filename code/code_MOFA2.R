#!/usr/bin/env Rscript
#This script executes MOFA2 and produces the corresponding file for the microbiome normalization
#the script takes 2 argument the data as produced by the simulate_data*.R scripts and the microbiome normalization

args = commandArgs(trailingOnly=TRUE)
library(MOFA2)

data = readRDS(args[1])
microbiome.norm = args[2]


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }



results.MOFA2 = lapply(1:1000, function(rep) {
	print(rep)

	species = data[[rep]]$Simulated.Microbiotes
        metabolites = data[[rep]]$Simulated.Metabolites
        
	species = data[[rep]]$Simulated.Microbiotes
        metabolites = data[[rep]]$Simulated.Metabolites

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

	Metabolites.rep= t(as.matrix(metabolites))
	Microbiome.rep = t(as.matrix(normalized.species))


	MOFA_obj = MOFA2::create_mofa(list("Metabo"=Metabolites.rep, "Micro"=Microbiome.rep))
	print("create mofa: done")
	MOFA_obj = MOFA2::prepare_mofa(MOFA_obj)
	print("prepare mofa: done")
	MOFA_obj = MOFA2::run_mofa(MOFA_obj, use_basilisk = TRUE)
	print("run mofa: done")
	MOFA2::get_variance_explained(MOFA_obj)

})

saveRDS(results.MOFA2, paste0("MOFA2_",microbiome.norm,"_Adenomas.RDS"))
