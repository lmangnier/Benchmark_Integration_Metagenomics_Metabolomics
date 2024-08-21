#!/usr/bin/env Rscript
#This script executes the Mantel method and produces the corresponding file for the microbiome normalization, correlation and distance kernel
#Different Distance kernels can be considered
#Metabolome data can be log transformed See the corresponding line in the code

#The script takes 5 arguments, 1) the microbiome normalization 2) the distance kernel 3) the correlation 4) data as produce by simulate_data*.R scripts 5) and the scenario either null or alternative


apply_nzv_filter <- function(x, freqCut = 95/5, uniqueCut = 10){
        nzv <- mixOmics::nearZeroVar(x = x, freqCut = freqCut, uniqueCut = uniqueCut)
        if(length(nzv$Position) > 0) x <- x[, -nzv$Position]
        return(x) }


args = commandArgs(trailingOnly=TRUE)

normalization = args[1]
dist.Metabolome = args[2]
method = args[3]

data = args[4]
scenario = args[5]

data = readRDS(data)
compute.Mantel = function(Microbiome,Metabolome, norm.Microbiome=c("clr","ilr", "alpha") ,dist.Metabolome=c("euclidean", "manhattan", "canberra"), method=c("spearman", "pearson")){
  
  
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
  
  #Compute euclidean distance on normalized microbiome data
  dist.microbiome = dist(normalized.microbiome, method = "euclidean")
  
  
  if(dist.Metabolome=="euclidean") dist.metabolome = dist(Metabolome, method = "euclidean")
  else if(dist.Metabolome=="canberra") dist.metabolome = dist(Metabolome, method = "canberra")
  else if(dist.Metabolome =="manhattan") dist.metabolome = dist(Metabolome, method = "manhattan")
  
  if(method=="pearson")  m = vegan::mantel(dist.microbiome, dist.metabolome, method="pearson")
  else if(method=="spearman")  m = vegan::mantel(dist.microbiome, dist.metabolome, method="spearman")
  
  return(m$signif)
}


pvalues = sapply(1:length(data), function(rep){
 print(rep)
 set.seed(runif(1,0,10000))
 species = data[[rep]]$Simulated.Microbiotes
 metabolites = data[[rep]]$Simulated.Metabolites #can be changed considering the log of the metabolites
 
 species = apply_nzv_filter(species)
 metabolites = apply_nzv_filter(metabolites)

 compute.Mantel(species, metabolites, norm.Microbiome=normalization,dist.Metabolome=dist.Metabolome, method=method)
})

saveRDS(pvalues, paste0("pvalues_Mantel_", scenario, "_", normalization, "_", dist.Metabolome,"_",method, "_filter.RDS")) 


