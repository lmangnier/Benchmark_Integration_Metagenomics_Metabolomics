#!/usr/bin/env Rscript
#This script executes the clr-lm model returning p-values 
#The script takes data as returned by the code_simulation*.R

args = commandArgs(trailingOnly=TRUE)

set.seed(1234)


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }

data = readRDS(args[1])

lm.results = lapply(1:length(data), function(x) {
        print(x)
        species = as.matrix(data[[x]]$Simulated.Microbiotes)
        metabolites = as.matrix(data[[x]]$Simulated.Metabolites) #log-transformation should be considered

        colnames(species) = paste0("Micro_", 1:ncol(species))
        colnames(metabolites) = paste0("Metabo_", 1:ncol(metabolites))

        species = apply_nzv_filter(species)
        metabolites = apply_nzv_filter(metabolites)

        clr.species = compositions::clr(species)

        list.results = list()

        for(col in colnames(metabolites)){
                list.results[[col]] = sapply(colnames(clr.species), function(x) {
                l = summary(lm(clr.species[,x]~metabolites[,col]))
                l$coefficients[2,4]
        })
        }


        list.results

})


saveRDS(lm.results,"results_lm_log.RDS")
