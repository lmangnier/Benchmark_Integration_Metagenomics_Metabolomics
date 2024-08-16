#This function simulates data based under the null for our Autism Scenario

#scenario == Autism: list.parameters should have 6 elements
#2 sparse correlation matrices
#lambda for Poisson 
#and the three parameters for the ZINB (mu, size and pstr0)

#list.parameters.Autism can be found in the input\Autism directory

list.parameters.Autism = readRDS("input\\Autism\\list_parameters_Autism.RDS")



simulate.realistic.data.alter.Autism = function(list.parameters, nrep=1000, nindiv=c(100,100), random.seed=1234){
  
  #nindiv is a vector where the first element is the number of individual in species, the second in metabolites
  
  #Sparse correlation matrices estimated by Spiec-Easi
  Cov.Species.Autism = list.parameters$Cov_species_Autism
  Cov.Metabolites.Autism = list.parameters$Cov_metabolites_Autism
  
  #Mean estimated by MLE 
  lambda.metabolites.Autism = list.parameters$lambda_metabolites
  
  #ZINB parameters estimated by MLE
  size.species.Autism = list.parameters$size_species
  prop.zeros.species.Autism = list.parameters$propzeros_species
  mu.species.Autism = list.parameters$mu_species
  
  set.seed(random.seed) #For reproducibility
  
  alter.rep.Autism = lapply(1: nrep, function(rep){
    print(rep)
    
    M1 = MASS::mvrnorm(nindiv[1], rep(0,ncol(Cov.Species.Autism)), Cov.Species.Autism)
    M2 = MASS::mvrnorm(nindiv[2], rep(0,ncol(Cov.Metabolites.Autism)), Cov.Metabolites.Autism)
    
    
    prop.associated.metabolites = round(runif(1,0.01,0.1),2)
    prop.associated.species = round(runif(1,0.01,0.1),2)
    
    prop.impacted.metabolites = round(runif(1,0.01,0.1),2)
    prop.impacted.species = round(runif(1,0.01,0.1),2)
    
    n.impacted.metabolites = round(ncol(Cov.Metabolites.Autism)*prop.impacted.metabolites)
    n.impacted.species = round(ncol(Cov.Species.Autism)*prop.impacted.species)
    
    n.associated.metabolites = ceiling(ncol(Cov.Metabolites.Autism)*prop.associated.metabolites)
    n.associated.species = ceiling(ncol(Cov.Species.Autism)*prop.associated.species)
    
    index.associated.species = sample(1:ncol(Cov.Species.Autism), n.associated.species)
    index.associated.metabolites = sample(1:ncol(Cov.Metabolites.Autism), n.associated.metabolites )
    
    index.impacted.species = sample(1:ncol(Cov.Species.Autism), n.impacted.species )
    index.impacted.metabolites = sample(1:ncol(Cov.Metabolites.Autism),  n.impacted.metabolites)
    
    vb=lapply(index.impacted.metabolites, function(x) sample(index.associated.species, sample(1:n.associated.species,1)))
    wb = lapply(vb, function(x) rnorm(length(x), 0,0.2))
    vn=lapply(index.impacted.species, function(x) sample(index.associated.metabolites, sample(1:n.associated.metabolites,1)))
    wn = lapply(vn, function(x) rnorm(length(x), 0,0.2))
    
    
    simulated.microbiotes = matrix(VGAM::qzinegbin(p = pnorm(M1), size = size.species.Autism, munb = mu.species.Autism, pstr0 = prop.zeros.species.Autism), ncol=ncol(Cov.Species.Autism),nrow=nindiv[1])
    simulated.metabolites = matrix(qpois(p = pnorm(M2), lambda = lambda.metabolites.Autism), ncol=ncol(Cov.Metabolites.Autism),nrow=nindiv[2])
    
    generate.assos.metabolites = sapply(1:length(index.impacted.metabolites), function(x) {
      mu = exp(simulated.microbiotes[,vb[[x]],drop=FALSE]%*%wb[[x]])
      mu[mu>1e+7] = 1e+7
      rpois(nrow(simulated.metabolites),mu)
    })
    
    generate.assos.species = sapply(1:length(index.impacted.species), function(x) {
      mu = exp(simulated.metabolites[,vn[[x]],drop=FALSE]%*%wn[[x]])
      mu[is.infinite(mu)] = max(mu[!is.infinite(mu)])
      
      mu[mu>1e+7] = 1e+7
      rnbinom(nrow(simulated.microbiotes), mu = mu, size=size.species.Autism[index.impacted.species[x]])
    })
    
    simulated.metabolites[, index.impacted.metabolites] = generate.assos.metabolites
    simulated.microbiotes[, index.impacted.species] = generate.assos.species
    
    colnames(simulated.microbiotes) = paste0("Micro", 1:ncol(simulated.microbiotes))
    colnames(simulated.metabolites) = paste0("Metabo", 1:ncol(simulated.metabolites))
    
    list("Simulated.Microbiotes"=simulated.microbiotes,"Simulated.Metabolites"=simulated.metabolites, "index.associated.species"=index.associated.species, "index.associated.metabolites"=index.associated.metabolites, "index.impacted.species"=index.impacted.species, "index.impacted.metabolites"=index.impacted.metabolites, "associations.metabolites"=vb, "associations.species"=vn)
  })
  
  return(alter.rep.Autism)
}
