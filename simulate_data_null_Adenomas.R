#This function simulates data based under the null for our Adenomas Scenario

#scenario == Adenomas: list.parameters should have 6 elements
#2 sparse correlation matrices
#mean for Gaussian 
#and the three parameters for the ZINB (mu, size and pstr0)

#list.parameters.Adenomas can be found in the input\Adenomas directory

list.parameters.Adenomas = readRDS("C:\\Users\\loicm\\Documents\\Paper_Benchmark\\Benchmark_Integration_Metagenomics_Metabolomics\\input\\Adenomas\\list_parameters_Adenomas.RDS")


simulate.realistic.data.null.Adenomas = function(list.parameters, nrep=1000, nindiv=100, random.seed=1234){
  
  #Sparse correlation matrices estimated by Spiec-Easi
  Cov.species.Adenomas = list.parameters$Cov_species_Adenomas
  Cov.Metabolites.Adenomas = list.parameters$Cov_metabolites_Adenomas
  
  #Mean estimated by MLE 
  mean.model.metabolites.Adenomas = list.parameters$mean_metabolites
  
  #ZINB parameters estimated by MLE
  size.species.Adenomas = list.parameters$size_species
  prop.zeros.species.Adenomas = list.parameters$propzeros_species
  mu.species.Adenomas = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  null.rep.Adenomas = lapply(1: nrep, function(rep){
    print(rep)
    
    multi.norm = MASS::mvrnorm(nindiv, rep(0,ncol(Cov.species.Adenomas)), Cov.species.Adenomas)
    simulated.metabolites = MASS::mvrnorm(nindiv, mean.model.metabolites.Adenomas, Cov.Metabolites.Adenomas)
    
    simulated.microbiotes = matrix(VGAM::qzinegbin(p = pnorm(multi.norm), size = size.species.Adenomas, munb = mu.species.Adenomas  ,pstr0 =  prop.zeros.species.Adenomas), ncol=ncol(Cov.species.Adenomas),nrow=nindiv)
    
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites)
    
  })
  
  return(null.rep.Adenomas)
  
}
