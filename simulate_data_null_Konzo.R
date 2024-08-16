#This function simulates data based under the null for our Konzo Scenario

#scenario == Konzo: list.parameters should have 5 elements
#2 sparse correlation matrices
#mean for Poisson
#and the two parameters for the NB (mu, size)

#list.parameters.Konzo can be found in the input\Konzo directory

list.parameters.Konzo = readRDS("C:\\Users\\loicm\\Documents\\Paper_Benchmark\\Benchmark_Integration_Metagenomics_Metabolomics\\input\\Konzo\\list_parameters_Konzo.RDS")

list.parameters.Konzo$lambda_metabolites
simulate.realistic.data.null.Konzo = function(list.parameters, nrep=1000, nindiv=100, random.seed=1234){
  
  #Sparse correlation matrices estimated by Spiec-Easi
  Cov.species.Konzo = list.parameters$Cov_species_Konzo
  Cov.Metabolites.Konzo = list.parameters$Cov_metabolites_Konzo
  
  #Mean estimated by MLE 
  lambda.metabolites.Konzo = list.parameters$lambda_metabolites
  
  #NB parameters estimated by MLE
  size.species.Konzo = list.parameters$size_species
  mu.species.Konzo = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  null.rep.Konzo = lapply(1: nrep, function(rep){
    print(rep)
    
    multi.norm = MASS::mvrnorm(nindiv, rep(0,ncol(Cov.species.Konzo)), Cov.species.Konzo)
    multi.norm1 = MASS::mvrnorm(nindiv, rep(0,ncol(Cov.Metabolites.Konzo)), Cov.Metabolites.Konzo)
    
    simulated.microbiotes = matrix(qnbinom(p = pnorm(multi.norm), size = size.species.Konzo, mu = mu.species.Konzo), ncol=ncol(Cov.species.Konzo),nrow=nindiv)
    simulated.metabolites = matrix(qpois(p = pnorm(multi.norm1), lambda = lambda.metabolites.Konzo), ncol=ncol(Cov.Metabolites.Konzo),nrow=nindiv)
    
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites)
    
  })
  
  return(null.rep.Konzo)
  
}

