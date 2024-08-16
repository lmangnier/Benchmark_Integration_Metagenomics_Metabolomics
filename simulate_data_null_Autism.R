#This function simulates data based under the null for our Autism Scenario

#scenario == Autism: list.parameters should have 6 elements
#2 sparse correlation matrices
#lambda for Poisson 
#and the three parameters for the ZINB (mu, size and pstr0)

#list.parameters.Autism can be found in the input\Autism directory

list.parameters.Autism = readRDS("input\\Autism\\list_parameters_Autism.RDS")

simulate.realistic.data.null.Autism = function(list.parameters, nrep=1000, nindiv=100, random.seed=1234){
  
  #Sparse correlation matrices estimated by Spiec-Easi
  Cov.species.Autism = list.parameters$Cov_species_Autism
  Cov.Metabolites.Autism = list.parameters$Cov_metabolites_Autism
    
  #Mean estimated by MLE 
  Poisson.model.metabolites.Autism = list.parameters$lambda_metabolites
    
  #ZINB parameters estimated by MLE
  size.species.Autism = list.parameters$size_species
  prop.zeros.species.Autism = list.parameters$propzeros_species
  mu.species.Autism = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  null.rep.Autism = lapply(1: nrep, function(rep){
    print(rep)
      
    multi.norm = MASS::mvrnorm(nindiv, rep(0,ncol(Cov.species.Autism)), Cov.species.Autism)
    multi.norm1 = MASS::mvrnorm(nindiv, rep(0,ncol(Cov.Metabolites.Autism)), Cov.Metabolites.Autism)
      
    simulated.metabolites = matrix(qpois(p = pnorm(multi.norm1), lambda=Poisson.model.metabolites.Autism), nrow=nindiv, ncol=ncol(Cov.Metabolites.Autism))
    simulated.microbiotes = matrix(VGAM::qzinegbin(p = pnorm(multi.norm), size = size.species.Autism, munb = mu.species.Autism  ,pstr0 =  prop.zeros.species.Autism), ncol=ncol(Cov.species.Autism),nrow=nindiv)
    
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites)
      
  })
  
  return(null.rep.Autism)
  
}

  